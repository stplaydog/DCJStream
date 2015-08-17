package gasts;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Collections;

public class Adjacency implements Serializable {
	public String name;
	public int num_gene, num_chr;
	public int[] reg_adj;
	public int p00 = 0;
	public ArrayList<Integer> cap_adj;

	public Adjacency() {
	}

	public Adjacency(Adjacency adj) {
		num_gene = adj.num_gene;
		num_chr = adj.num_chr;
		p00 = adj.p00;
		reg_adj = new int[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++)
			reg_adj[i] = adj.reg_adj[i];
		cap_adj = new ArrayList<Integer>(adj.cap_adj);
		name = adj.name;
	}

	public Adjacency(int n_g, int n_c, String nm) {
		num_gene = n_g;
		num_chr = n_c;
		reg_adj = new int[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++)
			reg_adj[i] = Constant.NULL;
		cap_adj = new ArrayList<Integer>(2 * num_chr);
		name = nm;
	}

	public int[] chromosomeCount() {
		int num_chro[] = new int[] { 0, 0, 0 };
		final int Linear = 0, Circular = 1, Piece = 2;
		boolean remains[] = new boolean[2 * num_gene];
		int genes[] = new int[2 * num_gene];
		for (int i = 0; i < num_gene; i++) {
			genes[2 * i] = 2 * i + 1;
			genes[2 * i + 1] = 2 * i;
		}

		for (int i = 0; i < 2 * num_gene; i++)
			remains[i] = true;

		int e_l, e_r;

		Collections.sort(cap_adj);
		for (Integer v : cap_adj) {
			if (!remains[v])
				continue;
			e_r = v;
			while (true) {
				e_l = genes[e_r];
				remains[e_l] = remains[e_r] = false;
				// looking for new adjacency
				e_r = reg_adj[e_l];
				if (e_r == Constant.CAP0) {
					num_chro[Linear]++;
					break;
				} else if (e_r == Constant.NULL) {
					num_chro[Piece]++;
					break;
				}
			}
		}

		// looking for pieces
		// now the remainings are either ciruclar chromosomes or pieces
		// System.out.println();
		for (int v = 0; v < 2 * num_gene; v++) {
			if (reg_adj[v] != Constant.NULL || !remains[v])
				continue;
			e_r = v;
			while (true) {
				e_l = genes[e_r];
				remains[e_l] = remains[e_r] = false;
				// looking for new adjacency
				e_r = reg_adj[e_l];
				if (e_r == Constant.NULL) {
					num_chro[Piece]++;
					break;
				}
			}
		}

		// looking for circular chromosomes
		// System.out.println();
		int start;
		for (int v = 0; v < 2 * num_gene; v++) {
			if (!remains[v])
				continue;
			start = v;
			e_r = v;
			while (true) {
				e_l = genes[e_r];
				remains[e_l] = remains[e_r] = false;
				// looking for new adjacency
				e_r = reg_adj[e_l];
				if (e_r == start) {
					num_chro[Circular]++;
					break;
				}
			}
		}
		num_chr = num_chro[Linear];
		return num_chro;
	}

	public int DCJ_distance(Adjacency adj) {
		return DCJ_distance(this, adj);
		// if(num_gene!=adj.num_gene) {
		// System.out.println("In calculating DCJ distance, two genomes contain unequal amount of genes!");
		// System.exit(1);
		// }
		// int n_chr=Math.max(num_chr, adj.num_chr);
		//
		// // calculating paths
		// int path_e=0, path_o=0, cycle=0;
		//
		// path_o+=(n_chr-cap_adj.size()/2); // null chromosomes;
		// path_o+=(n_chr-adj.cap_adj.size()/2); // null chromosomes;
		//
		// boolean remains[]=new boolean[2*num_gene];
		// for(int i=0;i<2*num_gene;i++) remains[i]=true;
		//
		// int e_l, e_r, start;
		// for(Integer v:cap_adj) {
		// e_r=v;
		// remains[e_r]=false;
		//
		// while(true) {
		// e_l=adj.reg_adj[e_r];
		// if(e_l==Constant.CAP0) {path_e++;break;}
		// remains[e_l]=false;
		//
		// e_r=reg_adj[e_l];
		// if(e_r==Constant.CAP0) {path_o++;break;}
		// remains[e_r]=false;
		// }
		// }
		//
		// for(Integer v:adj.cap_adj) {
		// if(!remains[v]) continue;
		// e_r=v;
		// remains[e_r]=false;
		//
		// while(true) {
		// e_l=reg_adj[e_r];
		// if(e_l==Constant.CAP0) {
		// System.out.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!");
		// System.exit(1);
		// }
		// remains[e_l]=false;
		//
		// e_r=adj.reg_adj[e_l];
		// if(e_r==Constant.CAP0) {path_o++;break;}
		// remains[e_r]=false;
		// }
		// }
		//
		// for(int v=0;v<2*num_gene;v++) {
		// if(!remains[v]) continue;
		// e_r=v;
		// start=v;
		// while(true) {
		// e_l=adj.reg_adj[e_r];
		// remains[e_l]=false;
		// remains[e_r]=false;
		// e_r=reg_adj[e_l];
		// if(e_r==start) {cycle++;break;}
		// }
		// }
		//
		// return num_gene+n_chr-cycle-path_e-path_o/2;
	}

	public static int DCJ_distance(Adjacency adj1, Adjacency adj2) {
		if (adj1.num_gene != adj2.num_gene) {
			System.out
					.println("In calculating DCJ distance, two genomes contain unequal amount of genes!");
			System.exit(1);
		}
		int n_gene = adj1.num_gene;
		int n_chr = Math.max(adj1.num_chr, adj2.num_chr);

		// if(adj1.num_chr!=adj1.cap_adj.size()/2 ||
		// adj2.num_chr!=adj2.cap_adj.size()/2)
		// System.out.println("bad in Adjacency in caclulating DCJ distance!");

		// calculating paths
		int path_e = 0, path_o = 0, cycle = 0;

		path_o += (n_chr - adj1.cap_adj.size() / 2); // null chromosomes;
		path_o += (n_chr - adj2.cap_adj.size() / 2); // null chromosomes;

		boolean remains[] = new boolean[2 * n_gene];
		for (int i = 0; i < 2 * n_gene; i++)
			remains[i] = true;

		int e_l, e_r, start;
		for (Integer v : adj1.cap_adj) {
			if (!remains[v])
				continue;
			e_r = v;
			remains[e_r] = false;

			while (true) {
				e_l = adj2.reg_adj[e_r];
				if (e_l == Constant.CAP0) {
					path_e++;
					break;
				}
				remains[e_l] = false;

				e_r = adj1.reg_adj[e_l];
				if (e_r == Constant.CAP0) {
					path_o++;
					break;
				}
				remains[e_r] = false;
			}
		}

		for (Integer v : adj2.cap_adj) {
			if (!remains[v])
				continue;
			e_r = v;
			remains[e_r] = false;

			while (true) {
				e_l = adj1.reg_adj[e_r];
				if (e_l == Constant.CAP0) {
					System.out
							.println("wrong in calculating DCJ distance, no cap edges in 'this' should remain!");
					System.exit(1);
				}
				remains[e_l] = false;

				e_r = adj2.reg_adj[e_l];
				if (e_r == Constant.CAP0) {
					path_o++;
					break;
				}
				remains[e_r] = false;
			}
		}

		for (int v = 0; v < 2 * n_gene; v++) {
			if (!remains[v])
				continue;
			e_r = v;
			start = v;
			while (true) {
				e_l = adj2.reg_adj[e_r];
				remains[e_l] = false;
				remains[e_r] = false;
				e_r = adj1.reg_adj[e_l];
				if (e_r == start) {
					cycle++;
					break;
				}
			}
		}

		int dcj = n_gene + n_chr - cycle - path_e - path_o / 2;
		return dcj;
	}

	public void findCommon(Adjacency adj) {
		if (num_gene != adj.num_gene) {
			System.out.println("different number of genes in findCommon!");
			System.exit(1);
		}
		name = "common_part";
		for (int i = 0; i < 2 * num_gene; i++) {
			if (reg_adj[i] != adj.reg_adj[i])
				reg_adj[i] = Constant.NULL;
		}

		cap_adj.retainAll(adj.cap_adj);

	}

	public Genome toGenome() {
		Genome gn = new Genome();
		gn.name = name;
		gn.num_chr = 0;
		num_chr = 0;
		gn.all_chromosomes = new ArrayList<String>();

		boolean remains[] = new boolean[2 * num_gene];
		int genes[] = new int[2 * num_gene];
		for (int i = 0; i < num_gene; i++) {
			genes[2 * i] = 2 * i + 1;
			genes[2 * i + 1] = 2 * i;
		}

		for (int i = 0; i < 2 * num_gene; i++)
			remains[i] = true;

		int e_l, e_r;

		Collections.sort(cap_adj);
		for (Integer v : cap_adj) {
			if (!remains[v])
				continue;
			e_r = v;
			num_chr++;
			StringBuilder ch = new StringBuilder();
			while (true) {
				e_l = genes[e_r];
				// e_r is the left side of the gene
				// e_l is the right side of the gene
				remains[e_l] = remains[e_r] = false;
				ch.append(label_to_gene(e_r) + " ");
				// looking for new adjacency
				e_r = reg_adj[e_l];
				if (e_r == Constant.CAP0) {
					break;
				} else if (e_r == Constant.NULL) {
					break;
				}
			}
			ch.append("$");
			gn.all_chromosomes.add(ch.toString());
		}

		// looking for pieces
		// now the remainings are either ciruclar chromosomes or pieces
		// System.out.println();
		for (int v = 0; v < 2 * num_gene; v++) {
			if (reg_adj[v] != Constant.NULL || !remains[v])
				continue;
			e_r = v;
			num_chr++;
			StringBuilder ch = new StringBuilder();
			while (true) {
				e_l = genes[e_r];
				remains[e_l] = remains[e_r] = false;
				ch.append(label_to_gene(e_r) + " ");
				// looking for new adjacency
				e_r = reg_adj[e_l];
				if (e_r == Constant.NULL) {
					break;
				}
			}
			ch.append("$");
			gn.all_chromosomes.add(ch.toString());
		}

		// looking for circular chromosomes
		// System.out.println();
		int start;
		for (int v = 0; v < 2 * num_gene; v++) {
			if (!remains[v])
				continue;
			start = v;
			e_r = v;
			num_chr++;
			StringBuilder ch = new StringBuilder();
			while (true) {
				e_l = genes[e_r];
				remains[e_l] = remains[e_r] = false;
				ch.append(label_to_gene(e_r) + " ");
				// looking for new adjacency
				e_r = reg_adj[e_l];
				if (e_r == start) {
					break;
				}
			}
			ch.append("$");
			gn.all_chromosomes.add(ch.toString());
		}

		gn.num_chr = num_chr;
		return gn;
	}

	int label_to_gene(int label) {
		// label is for the left end of the gene
		if (label % 2 == 0)
			return label / 2 + 1;
		else
			return -(label / 2 + 1);
	}
}
