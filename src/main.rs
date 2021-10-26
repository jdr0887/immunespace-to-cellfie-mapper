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

    let agent: ureq::Agent = ureq::AgentBuilder::new()
        .timeout_read(Duration::from_secs(10))
        .timeout_write(Duration::from_secs(10))
        .build();

    let gene_names = get_gene_names(&options.input)?;
    let gene_names = gene_names.into_iter().filter(|a| !a.is_empty()).collect_vec();
    let gene_symbol_name_to_hgnc_id_map = get_gene_name_to_hgnc_id_map(&agent, gene_names)?;

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

fn get_gene_name_to_hgnc_id_map(agent: &ureq::Agent, gene_names: Vec<String>) -> Result<HashMap<String, String>, Box<dyn error::Error>> {
    let mut map: HashMap<String, String> = HashMap::new();

    for chunk in gene_names.chunks(250) {
        let query = chunk.join("+OR+");
        let url = format!("http://rest.genenames.org/search/symbol/{}", query);
        let response = agent.get(url.as_str()).set("Accept", "application/json").call()?.into_string()?;
        let v: serde_json::Value = serde_json::from_str(response.as_str())?;
        v["response"]["docs"]
            .as_array()
            .unwrap()
            .into_iter()
            .map(|a| {
                let e: Entry = serde_json::from_value(a.clone()).unwrap();
                e
            })
            .for_each(|e| {
                map.entry(e.symbol).or_insert(e.hgnc_id.replace("HGNC:", ""));
            });
    }

    Ok(map)
}

fn get_gene_names(input: &path::PathBuf) -> Result<Vec<String>, Box<dyn error::Error>> {
    let file = fs::File::open(input)?;
    let reader = io::BufReader::new(file);
    let mut rdr = csv::Reader::from_reader(reader);

    let mut gene_names = vec![];
    for result in rdr.records().into_iter().skip(1) {
        let record = result?;
        let gene_name = record[0].to_string();
        gene_names.push(gene_name.replace("\"", ""));
    }
    gene_names.sort();
    gene_names.dedup();

    Ok(gene_names)
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

#[derive(Serialize, Deserialize, Debug)]
struct Entry {
    hgnc_id: String,
    score: f32,
    symbol: String,
}
