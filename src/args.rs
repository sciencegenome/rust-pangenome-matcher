use clap::Parser;

#[derive(Debug, Parser)]
#[clap(version)]

pub struct GenomeArgs {
    /// please provide the path to the first alignment file
    pub pafalignment1: String,
    /// please provide the path to the second alignment file
    pub pafalignment2: String,
}
