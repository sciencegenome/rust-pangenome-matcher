mod args;
use args::GenomeArgs;
use clap::Parser;
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
*Author Gaurav Sablok
*Universitat Potsdam
*Date 2024-12-11

rust-paftools: Building the entire paf tools in the rust for pangenome
construction from the metagenome and the pangenome and the pantranscriptomics.
This part is the pangenome matcher, and it takes two pangenome alignment, if you have aligned the
same query genome to the difference species then it will extract the regions aligned of the query
across the different species and will give you a tabulated information and also a harmonic mean
that how much pangenome works across the same species with the less divergence.

* */

fn main() {
    let args = GenomeArgs::parse();
    let paf_alignment = paf_alignments(&args.pafalignment1, &args.pafalignment2).unwrap();
     println!(
        "The pangenome matcher has produced the results: {}",
        paf_alignment
    );
}

fn paf_alignments(pathpaf1: &str, pathpaf2: &str) -> Result<String, Box<dyn Error>> {
    #[derive(Debug, Clone, PartialEq, PartialOrd)]
    struct QuerypafView {
        query: String,
        length: usize,
        start: usize,
        end: usize,
        strand: String,
    }

    #[derive(Debug, Clone, PartialEq, PartialOrd)]
    struct ReferencepafView {
        reference: String,
        length: usize,
        start: usize,
        end: usize,
        residuematch: usize,
        alignmentblock: usize,
    }

    let pafaligner1 = File::open(pathpaf1).expect("file not present");
    let pafreader1 = BufReader::new(pafaligner1);
    let mut pafquery1: Vec<QuerypafView> = Vec::new();
    let mut pafreference1: Vec<ReferencepafView> = Vec::new();
    for i in pafreader1.lines() {
        let line = i.expect("line not present");
        let linecheck = line.split("\t").collect::<Vec<_>>();
        pafquery1.push(QuerypafView {
            query: linecheck[0].to_string(),
                       length: linecheck[1].parse::<usize>().unwrap(),
                       start: linecheck[2].parse::<usize>().unwrap(),
                       end: linecheck[3].parse::<usize>().unwrap(),
                       strand: linecheck[4].to_string(),
        });
        pafreference1.push(ReferencepafView {
            reference: linecheck[5].to_string(),
                           length: linecheck[6].parse::<usize>().unwrap(),
                           start: linecheck[7].parse::<usize>().unwrap(),
                           end: linecheck[8].parse::<usize>().unwrap(),
                           residuematch: linecheck[9].parse::<usize>().unwrap(),
                           alignmentblock: linecheck[10].parse::<usize>().unwrap(),
        });
    }

    let pafaligner2 = File::open(pathpaf2).expect("file not present");
    let pafreader2 = BufReader::new(pafaligner2);
    let mut pafquery2: Vec<QuerypafView> = Vec::new();
    let mut pafreference2: Vec<ReferencepafView> = Vec::new();
    for i in pafreader2.lines() {
        let line = i.expect("line not present");
        let linecheck = line.split("\t").collect::<Vec<_>>();
        pafquery2.push(QuerypafView {
            query: linecheck[0].to_string(),
                       length: linecheck[1].parse::<usize>().unwrap(),
                       start: linecheck[2].parse::<usize>().unwrap(),
                       end: linecheck[3].parse::<usize>().unwrap(),
                       strand: linecheck[4].to_string(),
        });
        pafreference2.push(ReferencepafView {
            reference: linecheck[5].to_string(),
                           length: linecheck[6].parse::<usize>().unwrap(),
                           start: linecheck[7].parse::<usize>().unwrap(),
                           end: linecheck[8].parse::<usize>().unwrap(),
                           residuematch: linecheck[9].parse::<usize>().unwrap(),
                           alignmentblock: linecheck[10].parse::<usize>().unwrap(),
        });
    }

    #[derive(Debug, Clone, PartialEq, PartialOrd)]
    struct QueryComparativeHolder {
        query1: String,
        length1: usize,
        start1: usize,
        end1: usize,
        strand1: String,
        query2: String,
        length2: usize,
        start2: usize,
        end2: usize,
        strand2: String,
    }

    #[derive(Debug, Clone, PartialEq, PartialOrd)]
    struct RefComparativeHolder {
        reference1: String,
        length1: usize,
        start1: usize,
        end1: usize,
        reference2: String,
        length2: usize,
        start2: usize,
        end2: usize,
    }

    #[derive(Debug, Clone, PartialEq, PartialOrd)]
    struct ResidueAlignment {
        residuematch1: usize,
        residuematch2: usize,
        alignmentblock1: usize,
        alignmentblock2: usize,
    }

    let mut comparative_query: Vec<QueryComparativeHolder> = Vec::new();
    let mut comparative_reference: Vec<RefComparativeHolder> = Vec::new();
    let mut comparative_residue_alignment: Vec<ResidueAlignment> = Vec::new();
    for i in pafquery1.iter() {
        for j in pafquery2.iter() {
            for f1 in pafreference1.iter() {
                for f2 in pafreference2.iter() {
                    if i.query.to_string() == j.query.to_string() {
                        comparative_query.push(QueryComparativeHolder {
                            query1: i.query.clone(),
                                               length1: i.length,
                                               start1: i.start,
                                               end1: i.end,
                                               strand1: i.strand.clone(),
                                               query2: j.query.clone(),
                                               length2: j.length,
                                               start2: j.start,
                                               end2: j.end,
                                               strand2: j.strand.clone(),
                        });
                        comparative_reference.push(RefComparativeHolder {
                            reference1: f1.reference.clone(),
                                                   length1: f1.length,
                                                   start1: f1.start,
                                                   end1: f1.end,
                                                   reference2: f2.reference.clone(),
                                                   length2: f2.length,
                                                   start2: f2.start,
                                                   end2: f2.end,
                        });
                        comparative_residue_alignment.push(ResidueAlignment {
                            residuematch1: f1.residuematch,
                            residuematch2: f2.residuematch,
                            alignmentblock1: f1.alignmentblock,
                            alignmentblock2: f2.alignmentblock,
                        });
                    }
                }
            }
        }
    }

    let mut comparative_query_write =
    File::create("comparative_query.txt").expect("file not present");
    let mut comparative_ref_write =
    File::create("comparative_ref_write.txt").expect("file not present");

    println!("The comparative reference files have been written");
    for i in comparative_reference {
        writeln!(
            comparative_ref_write,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i.reference1, i.reference2, i.start1, i.start2, i.end1, i.end2, i.length1, i.length2
        )
        .expect("line not present");
    }

    println!("The comparative query files have been written");
    for i in comparative_query {
        writeln!(
            comparative_query_write,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            i.query1,
            i.length1,
            i.start1,
            i.end1,
            i.strand1,
            i.query2,
            i.length2,
            i.start2,
            i.end2,
            i.strand2
        )
        .expect("line nor present");
    }

    Ok(
        "The pangenome alignments have been parsed and the comparative results have been written"
        .to_string(),
    )
   }
