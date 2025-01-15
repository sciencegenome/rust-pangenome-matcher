# rust-pangenome-tools
 - writing all pangenome tools in rust. 
 - This part compares the two pangenome alignment with the same query and a different reference and anchors the query across both. 
 - please see the last commit message and if it says compiled binary then it is completed or else still in development version.

 ```
 cargo build

 ```
 ```
 λ gauravsablok rust-pangenome-matcher → λ git main* → ./target/debug/rust-pangenome-matcher -h
 Usage: rust-pangenome-matcher <PAFALIGNMENT1> <PAFALIGNMENT2>
 
 Arguments:
  <PAFALIGNMENT1>  please provide the path to the first alignment file
  <PAFALIGNMENT2>  please provide the path to the second alignment file

 Options:
  -h, --help     Print help
  -V, --version  Print version

 ```
 - to run the binary
 ```
λ gauravsablok rust-pangenome-matcher → λ git main* → ./target/debug/rust-pangenome-matcher ./sample-files/test1.paf ./sample-files/test2.paf
The comparative reference files have been written
The comparative query files have been written
The pangenome matcher has produced the results: The pangenome alignments have been parsed and the comparative results have been written

 ```
 Gaurav Sablok
