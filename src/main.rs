#[macro_use]
extern crate log;

use csv::Reader;
use humantime::format_duration;
use itertools::Itertools;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::error;
use std::fs;
use std::io;
use std::io::{BufRead, Read, Write};
use std::path;
use std::time::{Duration, Instant};
use structopt::StructOpt;

#[derive(StructOpt, Debug)]
#[structopt(name = "immunespace-to-cellfie-mapper", about = "fixes immunespace output to cellfie input")]
struct Options {
    #[structopt(short = "i", long = "input", long_help = "input", required = true, parse(from_os_str))]
    input: path::PathBuf,

    #[structopt(short = "m", long = "model", long_help = "model name", default_value = "MT_recon_2_2_entrez.mat")]
    model: String,

    #[structopt(short = "o", long = "output", long_help = "output", required = true, parse(from_os_str))]
    output: path::PathBuf,
}
fn main() -> Result<(), Box<dyn error::Error>> {
    let start = Instant::now();
    env_logger::init();

    let options = Options::from_args();
    debug!("{:?}", options);

    let genefilter_list = filter_genes_by_model(&options.model)?;
    let gene_symbol_name_to_hgnc_id_map = get_gene_symbol_name_to_hgnc_id_map(&genefilter_list)?;

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
                // warn!("problematic gene, look up individually: {}", gene);
                // let (_found_gene, hgnc_id) = get_gene_name_to_hgnc_id(&agent, gene.as_str())?;
                // let line = line.replace(gene.as_str(), hgnc_id.as_str()).replace("\"", "").replace(",", "\t");
                // writer.write_all(format!("{}\n", line).as_bytes())?;
            }
        }
    }

    info!("Duration: {}", format_duration(start.elapsed()).to_string());
    Ok(())
}

#[derive(Serialize, Deserialize, Debug)]
struct GeneDataEntry {
    symbol: String,
    hgnc_id: String,
    entrez_id: String,
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

fn filter_genes_by_model(model: &String) -> Result<Vec<String>, Box<dyn error::Error>> {
    let mut genes_filter_list: Vec<String> = vec![];

    let model_filter = include_str!("data/model-filter.csv");
    let mut rdr = csv::Reader::from_reader(model_filter.as_bytes());
    let headers = rdr.headers()?;
    let header_idx = headers
        .iter()
        .position(|h| h == model.replace(".mat", ""))
        .expect(format!("failed to get header index for: {}", model).as_str());

    for result in rdr.records().into_iter() {
        let record = result?;
        let value = record.get(header_idx);
        match value {
            None => {}
            Some(v) => {
                if !v.is_empty() {
                    genes_filter_list.push(v.to_string());
                }
            }
        }
    }
    genes_filter_list.sort();
    genes_filter_list.dedup();
    Ok(genes_filter_list)
}

fn get_gene_symbol_name_to_hgnc_id_map(genefilter_list: &Vec<String>) -> Result<HashMap<String, String>, Box<dyn error::Error>> {
    let genename_data = include_str!("data/genename-data.csv");

    let mut rdr = csv::Reader::from_reader(genename_data.as_bytes());

    let mut results_to_retain = vec![];
    for result in rdr.deserialize() {
        let record: GeneDataEntry = result?;
        if genefilter_list.contains(&record.entrez_id) {
            results_to_retain.push(record.symbol);
        }
    }

    debug!("results_to_retain: {:?}", results_to_retain);

    let tuples = genename_data
        .lines()
        .map(|line| {
            let split: Vec<String> = line.split(",").map(|s| s.to_string()).collect();
            (split.get(0).unwrap().clone(), split.get(1).unwrap().clone())
        })
        .collect_vec();

    let mut gene_symbol_name_to_hgnc_id_map: HashMap<String, String> = tuples.into_iter().collect();
    debug!("gene_symbol_name_to_hgnc_id_map.len(): {}", gene_symbol_name_to_hgnc_id_map.len());
    gene_symbol_name_to_hgnc_id_map.retain(|k, v| results_to_retain.contains(k));
    debug!("gene_symbol_name_to_hgnc_id_map.len(): {}", gene_symbol_name_to_hgnc_id_map.len());

    Ok(gene_symbol_name_to_hgnc_id_map)
}
