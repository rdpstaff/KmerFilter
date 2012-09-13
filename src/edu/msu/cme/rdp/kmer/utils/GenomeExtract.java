/*
 * Copyright (C) 2012 Jordan Fish <fishjord at msu.edu>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
package edu.msu.cme.rdp.kmer.utils;

import edu.msu.cme.rdp.eutils.EUtilsSeqHolder;
import edu.msu.cme.rdp.readseq.utils.gbregion.RegionParser;
import edu.msu.cme.rdp.readseq.utils.gbregion.SingleSeqRegion;
import edu.msu.cme.rdp.readseq.writers.FastaWriter;
import gov.nih.nlm.ncbi.www.soap.eutils.EFetchSequenceServiceStub.GBSet;
import java.io.File;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Map;
import javax.xml.stream.XMLInputFactory;

/**
 *
 * @author fishjord
 */
public class GenomeExtract {

    public static void main(String[] args) throws Exception {
        if (args.length != 2) {
            System.err.println("USAGE: GenomeExtract <genome_xml> <out_dir>");
            System.exit(1);
        }

        File[] genomeXmlFiles = new File(args[0]).listFiles(new FilenameFilter() {

            public boolean accept(File file, String string) {
                return string.contains(".xml");
            }
        });

        File outDir = new File(args[1]);
        outDir.mkdir();

        XMLInputFactory fact = XMLInputFactory.newInstance();

        for (File genomeXmlFile : genomeXmlFiles) {
            File genomeDir = new File(outDir, genomeXmlFile.getName().replace(".xml", ""));
            
            if(genomeDir.exists()) {
                System.err.println("Genome directory " + genomeDir + " exists, skipping extraction");
                continue;
            }
            
            System.err.println("Processing genome xml file " + genomeXmlFile);
            EUtilsSeqHolder seq = new EUtilsSeqHolder(GBSet.Factory.parse(fact.createXMLStreamReader(new FileReader(genomeXmlFile))).getGBSetSequence()[0].getGBSeq());
            
            
            if (!genomeDir.mkdir()) {
                throw new IOException("Failed to make directory " + genomeDir);
            }

            PrintStream geneMappingOut = new PrintStream(new File(genomeDir, seq.getPrimaryAccession() + "_gene_mapping.txt"));
            FastaWriter seqOut = new FastaWriter(new File(genomeDir, seq.getPrimaryAccession() + ".fasta"));
            FastaWriter protOut = new FastaWriter(new File(genomeDir, seq.getPrimaryAccession() + "_genes.fasta"));
            FastaWriter nuclOut = new FastaWriter(new File(genomeDir, seq.getPrimaryAccession() + "_genes_nucl.fasta"));

            String genomeSeq = seq.getSeqString();
            for (Map<String, String> qualifiers : seq.getFeatures("CDS", null, null)) {
                String gene = qualifiers.get("gene");
                if (gene == null) {
                    gene = qualifiers.get("product");
                }

                geneMappingOut.println(seq.getPrimaryAccession() + "\t" + gene + "\t" + qualifiers.get("location") + "\t" + qualifiers.get("codon_start") + "\t" + qualifiers.get("transl_table") + "\t" + qualifiers.get("location").contains("complement") + "\t" + qualifiers.get("protein_id"));

                if (qualifiers.containsKey("translation")) {
                    String protId = qualifiers.get("gene") + "_" + qualifiers.get("location");
                    if (qualifiers.containsKey("protein_id")) {
                        protId = qualifiers.get("protein_id");
                    }

                    protOut.writeSeq(protId, "gene=" + qualifiers.get("gene") + ",codon_start=" + qualifiers.get("codon_start") + ",transl_table=" + qualifiers.get("transl_table"), qualifiers.get("translation"));
                    
                    try {
                        SingleSeqRegion codedBy = RegionParser.parse(qualifiers.get("location"));
                        nuclOut.writeSeq(protId, ",gene=" + qualifiers.get("gene") + ",codon_start=" + qualifiers.get("codon_start") + ",transl_table=" + qualifiers.get("transl_table"), codedBy.getSeqRegion(genomeSeq));
                    } catch(Exception e) {
                        System.err.println("Unable to parse region for protId " + protId + " of genome " + seq.getPrimaryAccession() + ": " + e.getMessage());
                    }
                }
            }

            seqOut.writeSeq(seq.getPrimaryAccession(), genomeSeq);

            seqOut.close();
            geneMappingOut.close();
            protOut.close();
            nuclOut.close();
        }
    }
}
