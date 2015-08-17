package gasts;


import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;

public class Linearalization {
	static boolean has_random=false;
	GraphXu graph;
	Adjacency leaf_adj[];
	int info[][];
	int num_gene;
	int num_plasmids;
	int search[];
	Adjacency median;
	int min_max_order; // the lastest order of adding edges causing the first plasmid
	boolean get_initial;
	Random rd;

	public Linearalization (GraphXu g, Adjacency a1, Adjacency a2, Adjacency a3) {
		graph=g;
		num_gene=graph.num_gene;
		median=graph.median_adj;
		leaf_adj=new Adjacency[]{a1, a2, a3};
		info=new int[4][2*num_gene];
		min_max_order=Constant.NULL;
		get_initial=false;
		rd=new Random();

		get_linear_median();

	}

	public void get_linear_median() {

		while(true) {
			checkChromosomes();
//					System.out.printf("there are %d plasmids\n",num_plasmids);

			if(num_plasmids<=0 || median.cap_adj.size()==0 && num_plasmids<=1) {
				median.chromosomeCount();
//				System.out.printf("End of linearalization. There are %4d linear chromosoems and DCJ distance is %8d\n", median.num_chr,
//						Adjacency.DCJ_distance(median, leaf_adj[0])+Adjacency.DCJ_distance(median, leaf_adj[1])+Adjacency.DCJ_distance(median, leaf_adj[2]));
				return; 
			}

			for(int index=0;index<3;index++) 
				check_each_cycle(index);

			if(!get_initial) {
				get_initial=true;
				for(int i=0;i<graph.steps.size();i+=2) {
					int left=graph.steps.get(i);
					int right=graph.steps.get(i+1);
					if(left!=Constant.CAP0 && graph.adding_order[left]==min_max_order ||
							right!=Constant.CAP0 && graph.adding_order[right]==min_max_order) {
						int start_index=i;
						int adjacencies_from_linear=0; // to make sure there are enough adjacencies from linear chromosomes 
						for(int j=start_index+2;j<graph.steps.size();j+=2) {
							if(graph.steps.get(j)>=Constant.CAP0) continue;
							if(info[3][graph.steps.get(j)]>0) adjacencies_from_linear++;
							if(adjacencies_from_linear>=3*graph.num_plasmid) break;
						}
						if(adjacencies_from_linear<=3*graph.num_plasmid) {
							int j=start_index-2;
							while(adjacencies_from_linear<3*graph.num_plasmid && j>=0) {
								if(graph.steps.get(j)>=Constant.CAP0) continue;
								if(info[3][graph.steps.get(j)]>0) adjacencies_from_linear++;
								j-=2;
							}
							if(j<0) j=0;
							start_index=j;
						}
//						System.out.println("in Linearalization "+i+" "+start_index);
						search=new int[graph.steps.size()-start_index];
						for(int j=start_index;j<graph.steps.size();j++) search[j-start_index]=graph.steps.get(j);
						break;
					}
				}
			}
			//		System.out.printf("largest order %8d\n", graph.order);
			//		System.out.printf("size of steps array %8d\n",graph.steps.size());
			//		System.out.printf("size of search array %8d\n",search.length);

			int best_gain=-Constant.CAP0;
			ArrayList<EdgePair> best_pairs=null;
			for(int i=0;i<search.length;i+=2) {
				int left=search[i],right=search[i+1];

				if(left!=Constant.CAP0 && info[3][left]<0 || right!=Constant.CAP0 && info[3][right]<0) {
					
					boolean try_null=false;

//					if(left==Constant.CAP0 && right==Constant.CAP0) 
//						System.out.printf("Null chromosomes\n");
//					else if(left==Constant.CAP0 && right!=Constant.CAP0) 
//						System.out.printf("%8d %8d\t%4d %4d\t%4d " +
//								"%4d\n",-99, right, -99, info[3][right], -99, graph.adding_order[right]);
//					else if(left!=Constant.CAP0 && right==Constant.CAP0) 
//						System.out.printf("%8d %8d\t%4d %4d\t%4d " +
//								"%4d\n",left, -99, info[3][left], -99, graph.adding_order[left], -99);
//					else
//						System.out.printf("%8d %8d\t%4d %4d\t%4d " +
//								"%4d\n",left, right, info[3][left], info[3][right], graph.adding_order[left], graph.adding_order[right]);

					for(int j=0;j<search.length;j+=2) {
						if(j==i) continue;
						int l=search[j], r=search[j+1];
						if(l==Constant.CAP0 && r==Constant.CAP0) {
							if(try_null) continue;
							else {
								try_null=true;
								int gain=0;
								for(int index=0;index<3;index++) {
									if(info[index][left]<Constant.HOMO1) gain--;
									else if(info[index][left]>=Constant.HOMO2) gain++;
								}
								if(gain>best_gain) {
									best_gain=gain;
									best_pairs=new ArrayList<EdgePair>();
									best_pairs.add(new EdgePair(right, left, l, r,i,j));
								} 
								else if(gain==best_gain) {
									best_pairs.add(new EdgePair(right, left, l, r,i,j));
								}
								continue;
							}
						}
						if(l!=Constant.CAP0 && info[3][l]==info[3][left] || r!=Constant.CAP0 && info[3][r]==info[3][left]) continue;

						if(l==Constant.CAP0) { l=r;r=Constant.CAP0;}

						int gain1=0, gain2=0;
						for(int index=0;index<3;index++) {							
							if(info[index][left]==info[index][l] && info[index][left]<Constant.HOMO1) gain1++;
							else if(info[index][left]==-info[index][l] && info[index][left]<Constant.HOMO1) gain2++;
							else if(info[index][left]<Constant.HOMO1 && info[index][l]<Constant.HOMO1) {gain1--;gain2--;}
							else if(info[index][left]>=Constant.HOMO1 && info[index][left]<Constant.HOMO2 &&
									info[index][l]>=Constant.HOMO1 && info[index][l]<Constant.HOMO2) ;
							else if(info[index][left]>=Constant.HOMO2  && info[index][l]>=Constant.HOMO2) ;
							else if(info[index][left]>=Constant.HOMO1 && info[index][left]<Constant.HOMO2 && info[index][l]>=Constant.HOMO2 ||
									info[index][l]>=Constant.HOMO1 && info[index][l]<Constant.HOMO2 && info[index][left]>=Constant.HOMO2) {gain1++;gain2++;}
							else if(info[index][left]<Constant.HOMO1 && info[index][l]>=Constant.HOMO1 ||
									info[index][l]<Constant.HOMO1 && info[index][left]>=Constant.HOMO1) {gain1--;gain2--;}
							else 
							{
								System.out.println("Wrong in Liearalization, get_linear_median\n");
								System.out.printf("Is left larger than HOMO1? %b Larger than HOMO2? %b\n", info[index][left]>=Constant.HOMO1, info[index][left]>=Constant.HOMO2);
								System.out.printf("Is l larger than HOMO1? %b Larger than HOMO2? %b\n", info[index][l]>=Constant.HOMO1, info[index][l]>=Constant.HOMO2);
								System.exit(1);}
						}
						if(gain1>best_gain) {
							best_gain=gain1;
							best_pairs=new ArrayList<EdgePair>();
							best_pairs.add(new EdgePair(left, right, l, r,i,j));
						} 
						else if(gain1==best_gain) {
							best_pairs.add(new EdgePair(left, right, l, r,i,j));
						}
						if(gain2>best_gain) {
							best_gain=gain2;
							best_pairs=new ArrayList<EdgePair>();
							best_pairs.add(new EdgePair(right, left, l, r,i,j));
						} 
						else if(gain2==best_gain) {
							best_pairs.add(new EdgePair(right, left, l, r,i,j));
						}
					}
//					System.out.println();
				}
			}
			EdgePair selected_edgepair;
			if(has_random) selected_edgepair=best_pairs.get(rd.nextInt(best_pairs.size()));
			else selected_edgepair=best_pairs.get(0);
			//		System.out.printf("Best gain is %d, the two edges are (%4d,%4d) and (%4d,%4d)\n",
			//				best_gain, selected_edgepair.e1_l, selected_edgepair.e1_r, selected_edgepair.e2_l, selected_edgepair.e2_r);

			search[selected_edgepair.index1]=selected_edgepair.e1_l;search[selected_edgepair.index1+1]=selected_edgepair.e2_r;
			search[selected_edgepair.index2]=selected_edgepair.e2_l;search[selected_edgepair.index2+1]=selected_edgepair.e1_r;
			median.reg_adj[selected_edgepair.e1_l]=selected_edgepair.e2_r;
			median.reg_adj[selected_edgepair.e1_r]=selected_edgepair.e2_l;
			if(selected_edgepair.e2_l==Constant.CAP0 && selected_edgepair.e2_r==Constant.CAP0) {
				median.cap_adj.add(selected_edgepair.e1_l);
				median.cap_adj.add(selected_edgepair.e1_r);
			}

			else if(selected_edgepair.e2_r!=Constant.CAP0) {
				median.reg_adj[selected_edgepair.e2_l]=selected_edgepair.e1_r;
				median.reg_adj[selected_edgepair.e2_r]=selected_edgepair.e1_l;
			}
			else if(selected_edgepair.e2_r==Constant.CAP0) {
				median.reg_adj[selected_edgepair.e2_l]=selected_edgepair.e1_r;
				median.cap_adj.remove((Integer)selected_edgepair.e2_l);
				median.cap_adj.add(selected_edgepair.e1_l);
			}
		}
	}

	public void checkChromosomes() {
		int num_linear=0, num_circular=0;
		boolean remains[]=new boolean[2*num_gene];
		int genes[]=new int[2*num_gene];
		for(int i=0;i<num_gene;i++) {
			genes[2*i]=2*i+1;
			genes[2*i+1]=2*i;
		}

		for(int i=0;i<2*num_gene;i++) remains[i]=true;

		int e_l, e_r;

		Collections.sort(median.cap_adj);		
		for(Integer v:median.cap_adj) {
			if(!remains[v]) continue;
			num_linear++;
			e_r=v;
			while(true) {
				e_l=genes[e_r];
				remains[e_l]=remains[e_r]=false;
				info[3][e_l]=info[3][e_r]=num_linear;
				// looking for new adjacency
				e_r=median.reg_adj[e_l];
				if(e_r==Constant.CAP0) break;
			}
		}
		
		for(int i=0;i<2*num_gene;i++) {
			if(remains[i] && median.reg_adj[i]>=Constant.CAP0) System.out.println("Wrong! there are cap edges unchecked "+i);
			if(remains[i] && !remains[genes[i]]) System.out.println("Wrong with lableing! "+i);
		}

		// looking for circular chromosomes
		//		System.out.println();
		int start;
		for(int v=0;v<2*num_gene;v++) {
			if(!remains[v]) continue;
			num_circular--;
			int max_order=0;
			start=v;
			e_r=v;
			while(true) {
				e_l=genes[e_r];
				remains[e_l]=remains[e_r]=false;
				info[3][e_l]=info[3][e_r]=num_circular;
				if(!get_initial) {
					if(graph.adding_order[e_l]>max_order) max_order=graph.adding_order[e_l];
					if(graph.adding_order[e_r]>max_order) max_order=graph.adding_order[e_r];
				}
				// looking for new adjacency
				e_r=median.reg_adj[e_l];
				if(e_r==start) break;
			}
			if(!get_initial && max_order<min_max_order) min_max_order=max_order;
		}
		num_plasmids=-num_circular;
	}

	public void check_each_cycle(int index) {
		// 1, -1 are used for cycles and even paths
		// 2, -2 , 3, -3 are used for odd paths
		Adjacency adj1=median, adj2=leaf_adj[index];

		if(adj1.num_gene!=adj2.num_gene) {
			System.out.println("In checking cycles in Linearalization, two genomes contain unequal amount of genes!");
			System.exit(1);
		}

		int label=0, homo1=Constant.HOMO1, homo2=Constant.HOMO2;

		int n_gene=adj1.num_gene;
		// calculating paths

		boolean remains[]=new boolean[2*n_gene];
		for(int i=0;i<2*n_gene;i++) remains[i]=true;

		int e_l, e_r, start;
		for(Integer v:adj1.cap_adj) {
			if(!remains[v]) continue;
			int tmp_label;
			e_r=v;
			remains[e_r]=false;

			ArrayList<Integer> visited=new ArrayList<Integer>(); // to store vertices visited by going through this path
			visited.add(e_r);
			while(true) {
				e_l=adj2.reg_adj[e_r];
				if(e_l==Constant.CAP0) {tmp_label=++label;break;} // an even path
				remains[e_l]=false;
				visited.add(e_l);

				e_r=adj1.reg_adj[e_l];
				if(e_r==Constant.CAP0) {tmp_label=++homo1;break;} // an odd path with two cap edges with color median
				remains[e_r]=false;
				visited.add(e_r);
			}

			for(int i=0;i<visited.size();i+=2) info[index][visited.get(i)]=-tmp_label;
			for(int i=1;i<visited.size();i+=2) info[index][visited.get(i)]=tmp_label;
		}

		for(Integer v:adj2.cap_adj) { 
			// the only thing here is odd paths with two cap edges with color index
			// so 3 is used here
			if(!remains[v]) continue;
			e_r=v;
			remains[e_r]=false;
			homo2++;
			info[index][e_r]=homo2;

			while(true) {
				e_l=adj1.reg_adj[e_r];
				if(e_l==Constant.CAP0) {
					System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!"); 
					System.exit(1);
				}
				remains[e_l]=false;
				info[index][e_r]=-homo2;

				e_r=adj2.reg_adj[e_l];
				if(e_r==Constant.CAP0) break;
				remains[e_r]=false;
				info[index][e_r]=homo2;
			}
		}

		for(int v=0;v<2*n_gene;v++) {
			if(!remains[v]) continue;
			label++;
			e_r=v;
			start=v;
			while(true) {
				e_l=adj2.reg_adj[e_r];
				remains[e_l]=false;
				remains[e_r]=false;
				info[index][e_r]=-label;
				info[index][e_l]=label;
				e_r=adj1.reg_adj[e_l];
				if(e_r==start) break;
			}
		}
	}

}

class EdgePair {
	int e1_l, e1_r, e2_l, e2_r;
	int index1, index2;

	public EdgePair(int i1, int i2, int j1, int j2, int ind1, int ind2) {
		e1_l=i1;
		e1_r=i2;
		e2_l=j1;
		e2_r=j2;
		index1=ind1;
		index2=ind2;
	}
}
