import java.util.ArrayList;

public class VGPData {

	private String genomeId;
	private String speciesName;
	
	private int numSubreadBams = 0;
	private int numScrapsBams = 0;
	private int num10xR1 = 0;
	
	ArrayList<String> bionanoReBnx = new ArrayList<String>();
	ArrayList<String> bionanoReCmap = new ArrayList<String>();
	ArrayList<String> hicVenders = new ArrayList<String>();
	ArrayList<String> assemblies = new ArrayList<String>();
	
	private boolean updateCount = false;
	private int tech_count = 0;
	private int internal_score = 0;
	
	public VGPData(String genomeId, String speciesName) {
		this.genomeId = genomeId;
		this.speciesName = speciesName;
	}
	
	public void addPacbio(String file) {
		if (file.contains("subread") && file.endsWith(".bam")) {
			numSubreadBams++;
		} else if (file.contains("scraps") && file.endsWith(".bam")) {
			numScrapsBams++;
		}
		updateCount = true;
	}
	
	public void add10xR1(String file) {
		if (file.contains("r1")) {
			num10xR1++;
		}
		updateCount = true;
	}
	
	public void addHiC(String vendor) {
		if (!hicVenders.contains(vendor)) {
			hicVenders.add(vendor);
		}
		updateCount = true;
	}
	
	public void addBionano(String file) {
		ArrayList<String> bionanoRe = null;
		if (file.contains("bnx")) {
			bionanoRe = bionanoReBnx;
		} else if (file.contains("cmap")) {
			bionanoRe = bionanoReCmap;
		}
		
		if (file.contains("bspqi") && !bionanoRe.contains("BspQI")) {
			bionanoRe.add("BspQI");
		} else if (file.contains("bsssi") && !bionanoRe.contains("BssSI")) {
			bionanoRe.add("BssSI");
		} else if ((file.contains("dle") || file.contains("dls")) && !bionanoRe.contains("DLE1")) {
			bionanoRe.add("DLE1");
		}
		
		updateCount = true;
	}
	
	public void printVGPData(boolean isMDstyle) {
		updateTechCount();
		
		if (isMDstyle) {
			System.out.println("| " + speciesName + "\t" +
					"| " + genomeId + "\t" +
					"| " + tech_count + "\t" +
					"| " + numSubreadBams + "\t" + 
					"| " + numScrapsBams + "\t" +
					"| " + num10xR1 + "\t" +
					"| " + list(bionanoReBnx) + "\t" +
					"| " + list(bionanoReCmap) + "\t" +
					"| " + list(hicVenders) + "\t" +
					"| " + list(assemblies) + " |");
		} else {
			System.out.println(speciesName + "\t" +
					genomeId + "\t" +
					tech_count + "\t" +
					numSubreadBams + "\t" + 
					numScrapsBams + "\t" +
					num10xR1 + "\t" +
					list(bionanoReBnx) + "\t" +
					list(bionanoReCmap) + "\t" +
					list(hicVenders) + "\t" +
					list(assemblies));
		}
	}
	
	private void updateTechCount() {
		if (updateCount) {
			tech_count=0;
			internal_score++;
			if (numSubreadBams + numScrapsBams > 0) {
				tech_count++;
				internal_score++;
			}
			if (num10xR1 > 0) {
				tech_count++;
				internal_score++; 
			}
			if (bionanoReBnx.size() + bionanoReCmap.size() > 0) {
				tech_count++;
				if (bionanoReBnx.size() > 0) {
					internal_score ++;
				}
				if (bionanoReCmap.size() > 0) {
					internal_score ++;
				}
			}
			if (hicVenders.size() > 0) {
				tech_count++;
				internal_score += hicVenders.size();
			}

			// Give priorities to genomes with at least 1 assembly done
			if (assemblies.size() > 0) {
				tech_count ++;
				internal_score += 10;
			}
	
			// Earn extra +4 if all 4 platforms are collected
			if (tech_count >= 4) {
				internal_score += 4;
			}

			// Earn extra +6 if all 4 platforms + 1 assembly is collected
				if (tech_count >= 5) {
				internal_score += 6;
			}
		}
		
		updateCount = false;
	}
	
	public int getTechCount() {
		updateTechCount();
		return tech_count;
	}
	
	public int getInternalScore() {
		updateTechCount();
		return internal_score;
	}
	public static void printHeader(boolean isMDstyle) {
		if (isMDstyle) {
			System.out.println("| genome_name\t"
					+ "| species_id\t"
					+ "| tech_count\t"
					+ "| pacbio_subreads\t"
					+ "| pacbio_scrubs\t"
					+ "| 10x\t"
					+ "| bionano_bnx\t"
					+ "| bionano_cmap\t"
					+ "| hic \t"
					+ "| assembly |");
			System.out.println("| :---------- | :---------- | :---------- | :---------- | :---------- | :----- | :----- | :----- | :----- | :----- |");
		} else {
			System.out.println("genome_name\tspecies_id\t"
					+ "tech_count\t"
					+ "pacbio_subreads\tpacbio_scrubs\t10x\t"
					+ "bionano_bnx\tbionano_cmap\t"
					+ "hic\t"
					+ "assembly");
		}
	}
	
	private String list(ArrayList<String> list) {
		String elements = "";
		for(String element : list) {
			elements += element + ",";
		}
		if (!elements.equals("")) {
			elements = elements.substring(0, elements.length() - 1);
		}
		if (elements.equals("")) {
			elements = "X";
		}
		return elements;
	}

	public void addAssembly(String file) {
		if(file.contains(".fasta")) {
			assemblies.add(file.substring(file.indexOf("_") + 1, file.indexOf(".fasta")));
		}
		
	}
}
