extern crate clap;
extern crate noodles_bam;
extern crate noodles_sam;
extern crate noodles_fasta;

use std::io::{Write};
use clap::{Arg, App};
use std::collections::HashMap;
use std::{fs::File, io::BufReader};
use noodles_sam as sam;
use noodles_fasta as fasta;
use string_builder::Builder;
use noodles_sam::record::cigar::op::Kind;

fn main() -> std::io::Result<()> {
    let matches = App::new("SamToFastq")
        .version("1.0")
        .author("Aaron M. <aaron.mckenna@dartmouth.edu>")
        .about("Convert a SAM/BAM to fasta file, given a reference")
        .arg(Arg::with_name("input")
            .short("i")
            .long("input")
            .value_name("FILE")
            .required(true)
            .help("An input SAM/BAM file containing reads already aligned to the reference")
            .takes_value(true))
        .arg(Arg::with_name("reference")
            .short("r")
            .long("ref")
            .value_name("FILE")
            .required(true)
            .help("The reference we've align our reads to, can have multiple contigs")
            .takes_value(true))
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .value_name("FILE")
            .required(true)
            .help("the FASTA alignment output file")
            .takes_value(true))
        .arg(Arg::with_name("full_reference")
            .short("f")
            .long("full_reference")
            .help("should we try to output the full alignment, i.e. the gaps in the read from the beginning of the contig to the end"))
        .get_matches();

    // pull out the command line arguments and setup our ref, input, and output
    // ----------------------------------------------------------------
    let input_file = matches.value_of("input").unwrap();
    let output_file = matches.value_of("output").unwrap();
    let reference_file = matches.value_of("reference").unwrap();
    let full_reference: bool = matches.is_present("full_reference");

    // make a mapping of FASTA names to reference sequences
    // ----------------------------------------------------------------
    let mut reference_sequences = HashMap::new();

    println!("Loading reference sequences from {}...", reference_file);

    let mut fasta_reader = File::open(reference_file)
        .map(BufReader::new)
        .map(fasta::Reader::new)?;

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = match File::create(&output_file) {
        Err(why) => panic!("couldn't create {}: {}", output_file, why),
        Ok(file) => file,
    };

    for fasta_record in fasta_reader.records() {
        let fasta = fasta_record?;
        reference_sequences.insert(fasta.reference_sequence_name().to_string(), fasta.sequence().to_vec());
    }

    println!("Loaded {} reference sequences...", reference_sequences.len());

    // Setup our SAM/BAM reader and start parsing sequences
    // ----------------------------------------------------------------
    println!("Processing reads from {}...", input_file);

    match input_file {
        _x if _x.ends_with(".bam") || _x.ends_with(".BAM") => {
            //let reader = File::open(_x).map(bam::Reader::new)?;
            panic!("Not supported yet: {} ", input_file)
        }
        _x if _x.ends_with(".sam") || _x.ends_with(".SAM") => {
            let mut reader = File::open(_x).map(BufReader::new).map(sam::Reader::new)?;
            reader.read_header()?;

            for result in reader.records() {
                match result.map(|read| {
                    let sequence = read.sequence();
                    let reference_name = read.reference_sequence_name();
                    let cigar = read.cigar();
                    let alignment_start = read.position();

                    if let Some(ref_name) = reference_name {
                        if reference_sequences.contains_key(&ref_name.to_string()) && sequence.to_string() != "*" {
                        match alignment_start {
                            Some(p) => {
                                let position: usize = i32::from(p) as usize;
                                let mut reference_builder = Builder::default();
                                let mut read_builder = Builder::default();

                                let reference_sequence = &reference_sequences[&ref_name.to_string()];

                                // pad the front if requested
                                if full_reference {
                                    reference_builder.append(&reference_sequence[0..position]);
                                    read_builder.append(gap_of_length(position));
                                }

                                let mut current_ref_pos = position - 1;
                                let mut current_read_pos = 0;
                                for op in cigar.iter() {
                                    match op.kind() {
                                        Kind::Match => {
                                            let length: usize = op.len() as usize;
                                            let ref_seq = String::from_utf8(reference_sequence[current_ref_pos..(current_ref_pos + length)].to_vec());
                                            reference_builder.append(ref_seq.unwrap());
                                            current_ref_pos += length;
                                            let read_seq = sequence[current_read_pos..(current_read_pos + length)].iter().map(|c| char::from(*c)).collect::<String>();
                                            read_builder.append(read_seq);
                                            current_read_pos += length;
                                        }
                                        Kind::Insertion => {//Kind::Insertion => {
                                            let length: usize = op.len() as usize;
                                            reference_builder.append(gap_of_length(length));
                                            let read_seq = sequence[current_read_pos..(current_read_pos + length)].iter().map(|c| char::from(*c)).collect::<String>();
                                            read_builder.append(read_seq);
                                            current_read_pos += length;
                                        }
                                        Kind::Deletion => {
                                            let length: usize = op.len() as usize;
                                            let ref_seq = String::from_utf8(reference_sequence[current_ref_pos..(current_ref_pos + length)].to_vec());

                                            reference_builder.append(ref_seq.unwrap());
                                            current_ref_pos += length;
                                            let read_seq = gap_of_length(length);
                                            read_builder.append(read_seq);
                                        }
                                        Kind::SoftClip => {
                                            let length: usize = op.len() as usize;
                                            current_read_pos += length;
                                        }
                                        Kind::HardClip => {
                                            // these bases are already gone
                                        }
                                        _ => {
                                            panic!("Unsupported {} ", op);
                                        }
                                    };
                                };
                                // pad the front if requested
                                if full_reference {
                                    let remaining_reference = reference_sequence.len() - current_ref_pos;
                                    reference_builder.append(&reference_sequence[current_ref_pos..reference_sequence.len()]);
                                    read_builder.append(gap_of_length(remaining_reference));
                                }

                                writeln!(&mut file, ">{}", &ref_name.to_string()).expect("failed to write to file");
                                writeln!(&mut file, "{}", &reference_builder.string().unwrap()).expect("failed to write to file");
                                writeln!(&mut file, ">{}",read.read_name().unwrap()).expect("failed to write to file");
                                writeln!(&mut file, "{}", read_builder.string().unwrap()).expect("failed to write to file");
                            }
                            None => {}
                        };
                    };
                };
            })  {
                    Ok(_p) => { }
                    Err(_e) => panic!("Failed on matching"),
                };
        }
    }
    _ => panic!("We expect a file ending with .sam or .bam as input, we saw: {} ", input_file)
};

    // we're ok! <-- classic rust
    Ok(())
}


fn gap_of_length(x: usize) -> String {
    match String::from_utf8((0..x).map(|_| '-' as u8).collect::<Vec<u8>>()) {
        Ok(p) => { p }
        Err(_e) => panic!("Cant make a gap of length {}", x),
    }
}

