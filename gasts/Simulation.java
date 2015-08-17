package gasts;


import java.io.*;
import java.util.Random;

public class Simulation {
	int length;
	Newick<SimuNode> nw;
	SimuNode root;
	PrintWriter pw_leaf;
	PrintWriter pw_all, pw_real_tree;
	boolean regard_edge_length=false; 
	boolean use_random=false;
	int total_reversal=0;
	// if this true, then the edge lengths specified by the tree are used, the number of reversals is equal to edge lengths times length (the diameter)
	// otherwise, length is used as each edge length;
	int num_gene, num_chr;
	static int ii;
	boolean has_transposition=false;
	double transposition_rate=0;
	
	public Simulation(String tree, int n_gene, int n_chr, int edge) throws Exception  {
		length=edge;
//		pw_leaf=new PrintWriter(new File("leaf_"+tree+"_"+n_gene+"_"+n_chr+"_"+edge));
//		pw_all=new PrintWriter(new File("all_"+tree+"_"+n_gene+"_"+n_chr+"_"+edge));
		pw_leaf=new PrintWriter(new File("leaf_"+tree+"."+ii));
		pw_all=new PrintWriter(new File("all_"+tree+"."+ii));
		pw_real_tree=new PrintWriter(new File("real_"+tree+"."+ii));
		nw=new Newick<SimuNode>(new File(tree));
		nw.readTree(SimuNode.class);
		root=nw.reroot();
		num_gene=n_gene;
		num_chr=n_chr;
		root.edge_length=0;

	}
	
	public Simulation(String tree, int n_gene, int n_chr, int edge, boolean consider_edge_length) throws Exception  {
		this(tree, n_gene, n_chr, edge);
		regard_edge_length=consider_edge_length;
	}
	
	public Simulation(String tree, int n_gene, int n_chr, int edge, boolean consider_edge_length, boolean has_random, double rate_trans) throws Exception  {
		this(tree, n_gene, n_chr, edge);
		regard_edge_length=consider_edge_length;
		use_random=has_random;
		if(rate_trans>0) {
			has_transposition=true;
			transposition_rate=rate_trans;
		}
	}
	
	public void get_simulated_genomes ()  throws Exception {
		root.genome=new SimulatedGenome(num_gene,num_chr);
		pw_leaf.print(root.genome.get_string(root.name));
		pw_all.print(root.genome.get_string(root.name));
		simulate_each(root);
		pw_real_tree.println(root.get_Newick(true));
		pw_real_tree.println("total reversal distance "+total_reversal);
		pw_real_tree.close();
		pw_all.close();
		pw_leaf.close();
	}
	
	public void simulate_each(SimuNode node) throws Exception {
		for(int i=1;i<node.links.size();i++) {
			SimuNode child=(SimuNode) node.links.get(i);
			child.genome=new SimulatedGenome(node.genome);
			int num_event;
			if(regard_edge_length) {
				if(use_random) num_event=Math.max(1, Poisson.nextPoisson(child.edge_length*length));
				else num_event=(int) (child.edge_length*length+0.5);
//				else num_event=(int) (child.edge_length*(1-transposition_rate/1.6));
//				else num_event=(int) (child.edge_length*1.01);
			}
			else {
				if(use_random) num_event=Math.max(1, Poisson.nextPoisson(length));
				else num_event=(int) length;
			}
			
			if(!has_transposition) child.genome.reverse(num_event);
			else {
				Random rndm = new Random();
				for(int ii=0;ii<num_event;ii++) {
					if(rndm.nextDouble()<transposition_rate) {
						child.genome.tranpose(1);
					}
					else child.genome.reverse(1);
				}
			}
			int reversal=ReversalDistance.reversalDistance(node, child);
			child.edge_length=reversal;
			total_reversal+=reversal;
			if(child.isLeaf) pw_leaf.print(child.genome.get_string(child.name));
			pw_all.print(child.genome.get_string(child.name));
			simulate_each(child);
		}
	}
	
	public static void main(String[] args) throws Exception {
		if(args.length<6) {
			System.out.println("Please specify the tree file, number of genes, number of linear chromosomes, edge length or diameter, whether to use the edges lengths in the tree, whether to introduce randomness, the rate of transpositions!");
			System.exit(1);
		}
//		Simulation sm=null;
//		if(args.length==4)  sm=new Simulation(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]),Integer.parseInt(args[3]));
//		else if(args.length==5)  sm=new Simulation(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]),Integer.parseInt(args[3]), Boolean.parseBoolean(args[4]));
//		else if(args.length==6)  sm=new Simulation(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]),Integer.parseInt(args[3]), Boolean.parseBoolean(args[4]), Boolean.parseBoolean(args[5]));
//		sm.get_simulated_genomes();
		
		double rate_trans=0;
		if(args.length==7) rate_trans=Double.parseDouble(args[6]);
		
		for(int repeat=1;repeat<=10;repeat++) {
			ii=repeat;
			System.out.println("repeat="+ii);
			Simulation sm=null;
			do{
				sm=new Simulation(args[0], Integer.parseInt(args[1]), Integer.parseInt(args[2]), Integer.parseInt(args[3]),Boolean.parseBoolean(args[4]),Boolean.parseBoolean(args[5]),rate_trans);			
				sm.get_simulated_genomes();
				System.out.println(sm.total_reversal);
			}while(false);
//			while(sm.total_reversal<1068 || sm.total_reversal>1080);
//			}while(false);
		}
	}
}





