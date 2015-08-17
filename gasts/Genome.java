package gasts;


import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.ArrayList;
//this class should provide tools in reading file into Genome format or to generate random genomes
import java.util.Scanner;

public class Genome {
	String name;
	int num_chr; // assume there are only linear chromosomes
	ArrayList<String> all_chromosomes;

	public Genome() {}

	public static ArrayList<Genome> readFromFile(String file) throws Exception {
		BufferedReader br=new BufferedReader(new FileReader(new File(file)));
		String line;
		ArrayList<Genome> all_genomes=new ArrayList<Genome>();
		Genome gn=null;
		while((line=br.readLine())!=null) {
			//System.out.println(line);
			if(line.contains(">")) {
				if(gn!=null) {
					System.out.printf("Genome %s contains %4d chromosomes\n", gn.name, gn.num_chr);
					all_genomes.add(gn);
				}
				gn=new Genome();
				gn.name=line.replaceAll("[> ]", "") ;
				gn.num_chr=0;
				gn.all_chromosomes=new ArrayList<String>();
				continue;
			}
			else if(line.matches("\\s*?[-\\d].*?")) {
				//System.out.println(line);
				gn.all_chromosomes.add(line);
				gn.num_chr++;
				continue;
			}
		}
		if(gn!=null) {
			System.out.printf("Genome %s contains %4d chromosomes\n", gn.name, gn.num_chr);
			all_genomes.add(gn);
		}
		return all_genomes;
	}


	public Adjacency genome_to_adjacency (String linear_or_circular) {
		Adjacency adj=new Adjacency();
		adj.num_chr=num_chr;
		adj.name=name;

		ArrayList<ArrayList<Integer>> tmp_edges=new ArrayList<ArrayList<Integer>>(num_chr);
		// at this moment, we don't know number of genes
		// this block counts number of genes and store vertices into 2-dim array tmp_edges;
		adj.num_gene=0;
		if(linear_or_circular.equals("linear")) {
			for(String line:all_chromosomes) {
				ArrayList<Integer> tmp_chr=new ArrayList<Integer>();
				tmp_chr.add(Constant.CAP0);
				Scanner sc=new Scanner(line);
				while(sc.hasNextInt()) {

					adj.num_gene++;

					int gene=sc.nextInt();
					if(gene>0) {
						tmp_chr.add(2*gene-2);
						tmp_chr.add(2*gene-1);
					}
					else {
						tmp_chr.add(2*(-gene)-1);
						tmp_chr.add(2*(-gene)-2);
					}
				}	
				tmp_chr.add(Constant.CAP0);
				tmp_edges.add(tmp_chr);
			}
			
			System.out.println(name+"has "+adj.num_gene+" genes");

			// build adjacency arrays now
			adj.reg_adj=new int[2*adj.num_gene];
			adj.cap_adj=new ArrayList<Integer>(2*adj.num_chr);
			for(ArrayList<Integer> al: tmp_edges) {
				for(int index=0;index<al.size();index+=2) {
					int l=al.get(index), r=al.get(index+1);
					if(l!=Constant.CAP0 && r!=Constant.CAP0) {
						adj.reg_adj[l]=r;
						adj.reg_adj[r]=l;
					}
					else if(l==Constant.CAP0 && r==Constant.CAP0) adj.p00++;
					else if(l==Constant.CAP0 && r!=Constant.CAP0) {
						adj.reg_adj[r]=Constant.CAP0;
						adj.cap_adj.add(r);
					}
					else if(r==Constant.CAP0 && l!=Constant.CAP0) {
						adj.reg_adj[l]=Constant.CAP0;
						adj.cap_adj.add(l);
					}
				}
			}
		}
		else{
			for(String line:all_chromosomes) {
				ArrayList<Integer> tmp_chr=new ArrayList<Integer>();
				Scanner sc=new Scanner(line);
				if(!sc.hasNextInt()) continue;
				int first_gene=sc.nextInt();
				adj.num_gene++;
				int left, right;
				if(first_gene>0) {
					left=2*first_gene-2;
					right=2*first_gene-1;
				}
				else {
					left=-2*first_gene-1;
					right=-2*first_gene-2;
				}
				tmp_chr.add(right);
				while(sc.hasNextInt()) {

					adj.num_gene++;

					int gene=sc.nextInt();
					if(gene>0) {
						tmp_chr.add(2*gene-2);
						tmp_chr.add(2*gene-1);
					}
					else {
						tmp_chr.add(2*(-gene)-1);
						tmp_chr.add(2*(-gene)-2);
					}
				}	
				tmp_chr.add(left);
				tmp_edges.add(tmp_chr);
			}

			// build adjacency arrays now
			adj.reg_adj=new int[2*adj.num_gene];
			adj.cap_adj=new ArrayList<Integer>(2*adj.num_chr);
			for(ArrayList<Integer> al: tmp_edges) {
				for(int index=0;index<al.size();index+=2) {
					int l=al.get(index), r=al.get(index+1);
					adj.reg_adj[l]=r;
					adj.reg_adj[r]=l;
				}
			}
		}
		return adj;
	}	

	@Override public String toString() {
		StringBuilder sb=new StringBuilder();
		sb.append(">");sb.append(name);sb.append("\n");
		for(int i=0;i<all_chromosomes.size();i++) {
			sb.append("# chr");sb.append(i+1);sb.append("\n");
			sb.append(all_chromosomes.get(i));
			sb.append("\n");
		}
		sb.append("\n");

		return sb.toString();
	}
}
