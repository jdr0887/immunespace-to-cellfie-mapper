#[macro_use]
extern crate log;

use humantime::format_duration;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error;
use std::fs;
use std::io;
use std::io::{BufRead, Write};
use std::path;
use std::time::{Duration, Instant};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "immunespace-to-cellfie-mapper", about = "fixes immunespace output to cellfie input")]
struct Options {
    #[structopt(short = "i", long = "input", long_help = "input", required = true, parse(from_os_str))]
    input: path::PathBuf,

    #[structopt(short = "o", long = "output", long_help = "output", required = true, parse(from_os_str))]
    output: path::PathBuf,
}
fn main() -> Result<(), Box<dyn error::Error>> {
    let start = Instant::now();
    env_logger::init();

    let options = Options::from_args();
    debug!("{:?}", options);

    let genename_to_hgnc_id_data = include_str!("data/genename-to-hgnc-id.csv");

    let tuples = genename_to_hgnc_id_data
        .lines()
        .map(|line| {
            let split: Vec<String> = line.split(",").map(|s| s.to_string()).collect();
            (split.get(0).unwrap().clone(), split.get(1).unwrap().clone())
        })
        .collect_vec();

    let gene_symbol_name_to_hgnc_id_map: HashMap<String, String> = tuples.into_iter().collect();

    let agent: ureq::Agent = ureq::AgentBuilder::new()
        .timeout_read(Duration::from_secs(10))
        .timeout_write(Duration::from_secs(10))
        .build();

    let output_file = fs::File::create(options.output)?;
    let mut writer = io::BufWriter::new(output_file);

    let input_file = fs::File::open(options.input)?;
    let reader = io::BufReader::new(input_file);

    for line in reader.lines().skip(1) {
        let line = line?;
        let gene = line.split(",").next().unwrap().replace("\"", "");
        match gene_symbol_name_to_hgnc_id_map.get(gene.as_str()) {
            Some(gene_replacement) => {
                let line = line.replace(gene.as_str(), gene_replacement).replace("\"", "").replace(",", "\t");
                writer.write_all(format!("{}\n", line).as_bytes())?;
            }
            None => {
                warn!("problematic gene, look up individually: {}", gene);
                let (_found_gene, hgnc_id) = get_gene_name_to_hgnc_id(&agent, gene.as_str())?;
                let line = line.replace(gene.as_str(), hgnc_id.as_str()).replace("\"", "").replace(",", "\t");
                writer.write_all(format!("{}\n", line).as_bytes())?;
            }
        }
    }

    info!("Duration: {}", format_duration(start.elapsed()).to_string());
    Ok(())
}

#[derive(Serialize, Deserialize, Debug)]
struct Entry {
    hgnc_id: String,
    score: f32,
    symbol: String,
}

fn get_gene_name_to_hgnc_id(agent: &ureq::Agent, gene: &str) -> Result<(String, String), Box<dyn error::Error>> {
    let url = format!("http://rest.genenames.org/search/prev_symbol/{}", gene);
    let response = agent.get(url.as_str()).set("Accept", "application/json").call()?.into_string()?;
    let v: serde_json::Value = serde_json::from_str(response.as_str())?;
    let mut ret = vec![];
    v["response"]["docs"]
        .as_array()
        .unwrap()
        .into_iter()
        .map(|a| {
            let e: Entry = serde_json::from_value(a.clone()).unwrap();
            e
        })
        .for_each(|e| ret.push((e.symbol, e.hgnc_id.replace("HGNC:", ""))));
    Ok(ret.first().unwrap().clone())
}
