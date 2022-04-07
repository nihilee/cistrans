use std::collections::HashMap;
use std::env;
use std::error::Error;
use csv::StringRecord;
use itertools::Itertools;
use rust_htslib::{bam, bam::Read, bam::IndexedReader};
use std::str;

// static GENOME: &'static str = "/data/RefData/Genome/hg19/ucsc.hg19.fa";

fn main() {
    // println!("Hello, world!");

    let args: Vec<String> = env::args().collect::<Vec<_>>();

    let v: Vec<&str> = args.iter().map(|x| &**x).collect();

    // println!("{:?}", v);
    cis_tran(&v);

}

//
#[derive(Debug)]
struct Reads {
    ins: Vec<String>,
    del: Vec<String>,
    snp: Vec<String>,
    count: i32
}

fn haplotypes<'a>(site_f:&'a str, site_b:&'a str, bam:&'a str) -> Result<Vec<&'a str>, Box<dyn Error>> {
    let site_f:Vec<_> = site_f.split(":").collect();
    let site_b:Vec<_> = site_b.split(":").collect();

    let mut bam = IndexedReader::from_path(bam)?;
    bam.fetch((site_f[0], site_f[1].parse::<i64>()?, site_b[1].parse::<i64>()?));

    let mut site_f_reads = Reads {ins:Vec::new(), del:Vec::new(), snp:Vec::new(), count:0};
    let mut site_b_reads = Reads {ins:Vec::new(), del:Vec::new(), snp:Vec::new(), count:0};

    for p in bam.pileup(){
        let pileup = p.unwrap();

        if pileup.pos()+1 == site_f[1].parse::<u32>().unwrap() {
            // println!("pos is {} , {}", site_f[1], pileup.pos());
            for aln in pileup.alignments(){
                match aln.indel() {
                    bam::pileup::Indel::Ins(_) => {
                        // println!("ins: {}", len);
                        site_f_reads.ins.push(
                            str::from_utf8(aln.record().qname())?.to_string()
                        );
                        site_f_reads.count += 1;
                    },
                    bam::pileup::Indel::Del(_) => {
                        // println!("del: {}", len);
                        site_f_reads.del.push(
                            str::from_utf8(aln.record().qname())?.to_string()
                        );
                        site_f_reads.count += 1;
                    },
                    bam::pileup::Indel::None => {
                        if  !aln.is_del() && !aln.is_refskip() && str::from_utf8( &[aln.record().seq()[aln.qpos().unwrap()]] ).unwrap() == site_f[2].split("_").collect::<Vec<_>>()[1]{
                            site_f_reads.snp.push(
                                str::from_utf8(aln.record().qname())?.to_string()
                            );
                            site_f_reads.count += 1;
                        }
                    }
                }
            }
        } else if pileup.pos()+1 == site_b[1].parse::<u32>().unwrap() {
            // println!("pos is {} , {}", site_b[1], pileup.pos());
            for aln in pileup.alignments(){
                match aln.indel() {
                    bam::pileup::Indel::Ins(_) => {
                        // println!("ins: {}", len);
                        site_b_reads.ins.push(
                            str::from_utf8(aln.record().qname())?.to_string()
                        );
                        site_b_reads.count += 1;
                    },
                    bam::pileup::Indel::Del(_) => {
                        // println!("del: {}", len);
                        site_b_reads.del.push(
                            str::from_utf8(aln.record().qname())?.to_string()
                        );
                        site_b_reads.count += 1;
                    },
                    bam::pileup::Indel::None => {
                        if  !aln.is_del() && !aln.is_refskip() && str::from_utf8( &[aln.record().seq()[aln.qpos().unwrap()]] ).unwrap() == site_b[2].split("_").collect::<Vec<_>>()[1]{
                            site_b_reads.snp.push(
                                str::from_utf8(aln.record().qname())?.to_string()
                            );
                            site_b_reads.count += 1;
                        }
                    }
                }
            }
        }
    }

    let common = common_reads(&site_f_reads, &site_b_reads) as f64;
    let cis_f =  common / site_f_reads.count as f64;
    let trans_f = 1.0 - cis_f;
    let cis_b = common / site_b_reads.count as f64;
    let trans_b = 1.0 - cis_b;

    // println!("commom={}, site_f_reads={}, site_b_reads={}", common, site_f_reads.count, site_b_reads.count);
    // println!("cis_f={:.6},cis_b={:.6},trans_f={:.6},trans_b={:.6}", cis_f, cis_b, trans_f, trans_b);

    if cis_f < 0.15 && cis_b < 0.15 {
        Ok(vec!["trans"])
    } else if cis_f > 0.15 && cis_b > 0.15 && trans_f > 0.15 && trans_b > 0.15 {
        Ok(vec!["trans", "cis"])
    } else {
        Ok(vec!["cis"])
    }

}

fn common_reads(r1:&Reads, r2:&Reads) -> i32{
    let mut common = 0;
    let r1_del = r1.del.clone();
    let r1_ins = r1.ins.clone();
    let r1_snp = r1.snp.clone();

    let r2_del = r2.del.clone();
    let r2_ins = r2.ins.clone();
    let r2_snp = r2.snp.clone();

    for i in r1_del {
        if r2_del.contains(&i) {
            common += 1;
        }
    }

    for i in r1_ins {
        if r2_ins.contains(&i) {
            common += 1;
        }
    }

    for i in r1_snp {
        if r2_snp.contains(&i) {
            common += 1;
        }
    }

    common
}

fn add_cistrans(bam: &str, all: &mut HashMap<String, HashMap<String, String>>) -> Result<HashMap<String, HashMap<String, String>> ,Box<dyn Error>> {
    let mut pair = vec![];

    let all_shadow = all.to_owned();

    for item in all_shadow.iter().combinations(2) {
        let dist = (item[0].1.get("POS").unwrap().parse::<i32>().unwrap() - item[1].1.get("POS").unwrap().parse::<i32>().unwrap() ).abs();
        if item[0].1.get("CHR") == item[1].1.get("CHR") && dist > 0 && dist < 25 {
            pair.push((item[0].0, item[1].0));
        }
    }

    // println!("{:?}", pair);

    for (k, pk) in pair {
        let i = all.get(k).unwrap();
        let pi = all.get(pk).unwrap();

        let ct = haplotypes(pk, k, bam)?;

        let mut ct_type = vec![];
        let mut pre_ct_type = vec![];

        for t in ct {
            ct_type.push(format!("{}:{}({})", t, pi.get("HGVS.p").unwrap(), pi.get("HGVS.c").unwrap()));
            pre_ct_type.push(format!("{}:{}({})", t, i.get("HGVS.p").unwrap(), i.get("HGVS.c").unwrap()));
        }
        let ct_type = ct_type.join(";");
        let pre_ct_type = pre_ct_type.join(";");

        if all.get(k).unwrap().get("cistrans").unwrap() == "." {
            all.get_mut(k).unwrap().insert("cistrans".to_string(), ct_type);
        } else {
            let ty = all.get_mut(k).unwrap().entry("cistrans".to_string()).or_default();
            *ty = format!("{};{}", ty, ct_type);
        }

        if all.get(pk).unwrap().get("cistrans").unwrap() == "." {
            all.get_mut(pk).unwrap().insert("cistrans".to_string(), pre_ct_type);
        } else {
            let ty = all.get_mut(pk).unwrap().entry("cistrans".to_string()).or_default();
            *ty = format!("{};{}", ty, pre_ct_type);
        }
    }

    Ok(all.clone())
}


fn cis_tran(args: &[&str]) -> Result<(), Box<dyn Error>> {
    let bam = args[1];

    if args.len() == 3 {
        let f1 = args[2];
        let (h, mut a) = get_alldict(f1).unwrap();
        let new_a = add_cistrans(bam, &mut a)?;
        let mut wtr = csv::WriterBuilder::new().delimiter(b'\t').from_path(f1)?;

        wtr.write_record(&h)?;
        for rec in new_a.values(){
            let mut v = vec![];
            for i in h.iter() {
                v.push(rec.get(i).unwrap());
            }
            wtr.write_record(v)?;
        }
        wtr.flush()?;

    } else if args.len() == 4 {
        let f1 = args[2];
        let f2 = args[3];

        // println!("{}, {}", f1, f2);
        let (h1, a1) = get_alldict(f1).unwrap();
        let (h2, a2) = get_alldict(f2).unwrap();
        let mut all_dict = a1.clone().into_iter().chain(a2.clone()).collect();

        all_dict = add_cistrans(bam, &mut all_dict)?;

        let mut new_a1 = HashMap::new();
        let mut new_a2 = HashMap::new();

        for i in a1.keys() {
            new_a1.insert(i, all_dict.get(i).unwrap());
        }
        for i in a2.keys() {
            new_a2.insert(i, all_dict.get(i).unwrap());
        }

        let mut wtr1 = csv::WriterBuilder::new().delimiter(b'\t').from_path(f1)?;
        wtr1.write_record(&h1)?;
        for rec in new_a1.values(){
            let mut v = vec![];
            for i in h1.iter() {
                v.push(rec.get(i).unwrap());
            }
            wtr1.write_record(v)?;
        }
        wtr1.flush()?;

        let mut wtr2 = csv::WriterBuilder::new().delimiter(b'\t').from_path(f2)?;
        wtr2.write_record(&h2)?;
        for rec in new_a2.values(){
            let mut v = vec![];
            for i in h2.iter() {
                v.push(rec.get(i).unwrap());
            }
            wtr2.write_record(v)?;
        }
        wtr2.flush()?;
    }

    Ok(())
}

fn get_alldict(infile: &str) -> Result<(StringRecord, HashMap<String, HashMap<String, String>>), Box<dyn Error>> {
    let mut all: HashMap<String, HashMap<String, String>> = HashMap::new();

    let risk: HashMap<&str, u8> = HashMap::from(
        [("Pathogenic", 1),
            ("Likely pathogenic", 2),
            ("Uncertain significance", 3),
            ("Likely benign", 4),
            ("Benign", 5)]);

    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(true)
        .delimiter(b'\t')
        .from_path(infile)?;

    let header = rdr.headers()?.clone();

    for result in rdr.records() {
        let record = result?;

        let mut tmp: HashMap<_, _> =
            header.iter()
                .map(|s| s.to_string())
                .zip(record.iter().map(|s| s.to_string()))
                .collect();
        tmp.insert("cistrans".to_string(), ".".to_string());

        let k = format!("{}:{}:{}_{}",
                        tmp.get("CHR").unwrap(),
                        tmp.get("POS").unwrap(),
                        tmp.get("REF").unwrap(),
                        tmp.get("ALT").unwrap()
        );

        if all.contains_key(&k) {
            let a = risk.get("InterVar")
                .unwrap_or(&9);
            let b = risk.get(all.get(&k).unwrap().get("InterVar").unwrap().as_str())
                .unwrap_or(&9);
            if a < b {
                all.insert(k, tmp);
            }
        } else {
            all.insert(k, tmp);
        }
    }

    Ok((header, all.clone()))
}