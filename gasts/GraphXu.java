package gasts;

import java.util.ArrayList;
import java.io.*;
import java.util.*;

public class GraphXu implements Serializable {
	public int neighbor[][];
	ArrayList<ArrayList<Integer>> cap0_adj;
	ArrayList<ArrayList<Integer>> cap1_adj;
	int num_gene, num_chr, size, num_reg_vtx, num_free_caps;
	public int cycle_number;
	int path11[], path00[];
	int upper_bound;
	int lower_bound;
	int num_plasmid;

	int adding_order[]; // store the precedent order of adding edges. The
						// earliest added edges are thought most confident
	int order;
	boolean confident; // if AS0 or AS2 m=2 are meet, confident=false;

	boolean remains[];
	int new_label[], map[];
	Adjacency median_adj;

	ArrayList<Integer> steps; // to record the added edges in their precedence
								// order

	public void graph_vis_all() {
		System.out.print("graph G {\n");
		for (int i = 0; i < num_gene * 2; i++) {
			if (neighbor[i][0] > i)
				System.out.printf("%d -- %d [color=red];\n", i, neighbor[i][0]);
			if (neighbor[i][1] > i)
				System.out
						.printf("%d -- %d [color=blue];\n", i, neighbor[i][1]);
			if (neighbor[i][2] > i)
				System.out.printf("%d -- %d [color=green];\n", i,
						neighbor[i][2]);
		}
		System.out.print("}\n");
	}

	public GraphXu(Adjacency adj1, Adjacency adj2, Adjacency adj3) {
		Adjacency[] all_adj = new Adjacency[] { adj1, adj2, adj3 };
		// System.out.printf("number of chromosomes %4d %4d %4d\n",
		// adj1.num_chr,adj2.num_chr,adj3.num_chr);
		if (all_adj[0].num_gene != all_adj[1].num_gene
				|| all_adj[0].num_gene != all_adj[2].num_gene) {
			System.out
					.println("The three input genomes have different gene content!");
			System.exit(1);
		}
		num_gene = all_adj[0].num_gene;
		num_chr = Math.max(all_adj[0].num_chr,
				Math.max(all_adj[1].num_chr, all_adj[2].num_chr));
		size = num_gene + num_chr;
		num_reg_vtx = 2 * num_gene;
		num_free_caps = 2 * num_chr;

		remains = new boolean[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++)
			remains[i] = true;
		new_label = new int[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++)
			new_label[i] = i;
		map = new int[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++)
			map[i] = i;
		median_adj = new Adjacency(num_gene, num_chr, "median");

		neighbor = new int[num_reg_vtx][3];
		cap0_adj = new ArrayList<ArrayList<Integer>>(3);
		cap1_adj = new ArrayList<ArrayList<Integer>>(3);
		path00 = new int[3];
		path11 = new int[] { 0, 0, 0 };

		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < num_reg_vtx; j++)
				neighbor[j][i] = all_adj[i].reg_adj[j];
			cap0_adj.add(new ArrayList<Integer>(all_adj[i].cap_adj));
			cap1_adj.add(new ArrayList<Integer>());
			path00[i] = all_adj[i].p00 + num_chr - all_adj[i].num_chr;
		}

		adding_order = new int[2 * num_gene];
		order = 1;
		confident = true;

		steps = new ArrayList<Integer>();
		get_bounds();
	}

	public GraphXu(GraphXu g) {
		// copying numbers
		num_gene = g.num_gene;
		num_chr = g.num_chr;
		size = g.size;
		num_reg_vtx = g.num_reg_vtx;
		num_free_caps = g.num_free_caps;
		cycle_number = g.cycle_number;
		upper_bound = g.upper_bound;
		lower_bound = g.lower_bound;
		num_plasmid = g.num_plasmid;

		// 1-dim of size 3 arrays
		path00 = new int[3];
		path11 = new int[3];
		for (int c = 0; c < 3; c++) {
			path11[c] = g.path11[c];
			path00[c] = g.path00[c];
		}

		remains = new boolean[2 * num_gene];
		new_label = new int[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++) {
			remains[i] = g.remains[i];
			new_label[i] = g.new_label[i];
		}
		map = new int[num_reg_vtx];
		for (int i = 0; i < num_reg_vtx; i++)
			map[i] = g.map[i];
		median_adj = new Adjacency(g.median_adj);

		// copy adjacencies (edges)
		neighbor = new int[num_reg_vtx][3];
		for (int i = 0; i < num_reg_vtx; i++)
			for (int c = 0; c < 3; c++)
				neighbor[i][c] = g.neighbor[i][c];

		cap0_adj = new ArrayList<ArrayList<Integer>>(3);
		cap1_adj = new ArrayList<ArrayList<Integer>>(3);
		for (int c = 0; c < 3; c++) {
			cap0_adj.add(new ArrayList<Integer>(g.cap0_adj.get(c)));
			cap1_adj.add(new ArrayList<Integer>(g.cap1_adj.get(c)));
		}

		adding_order = new int[2 * num_gene];
		for (int i = 0; i < 2 * num_gene; i++)
			adding_order[i] = g.adding_order[i];
		order = g.order;
		confident = g.confident;

		steps = new ArrayList<Integer>(g.steps);
	}

	public void shrink(int left, int right) {
		// to shrink CMBG graph.
		// it can be used in constructing CMBG when the median genome is already
		// partially assigned
		// or to be used in the median solver to reduece grpah size

		// depending whether left, right take value of Constant.CAP0,
		// Constant.CAP1, NULL, there are a lot decisions here.

		size--;

		if (left != Constant.CAP0 && right != Constant.CAP0
				&& left != Constant.CAP1 && right != Constant.CAP1) {
			adding_order[map[left]] = adding_order[map[right]] = order;
			median_adj.reg_adj[map[left]] = map[right];
			median_adj.reg_adj[map[right]] = map[left];
			for (int c = 0; c < 3; c++) {
				int l_c = neighbor[left][c];
				int r_c = neighbor[right][c];
				// a lot things to consider, depending whether l_c, r_c take
				// values Constant.CAP0, Constant.CAP1, NULL

				if (l_c != Constant.CAP0 && l_c != Constant.CAP1
						&& l_c != Constant.NULL) {
					if (r_c != Constant.CAP0 && r_c != Constant.CAP1
							&& r_c != Constant.NULL) {
						if (r_c == left && l_c == right)
							cycle_number++;
						neighbor[r_c][c] = l_c;
						neighbor[l_c][c] = r_c;
					} else if (r_c == Constant.CAP0) {
						cap0_adj.get(c).remove((Integer) right);
						cap0_adj.get(c).add(l_c);
						neighbor[l_c][c] = Constant.CAP0;
					} else if (r_c == Constant.CAP1) {
						cap1_adj.get(c).remove((Integer) right);
						cap1_adj.get(c).add(l_c);
						neighbor[l_c][c] = Constant.CAP1;
					} else if (r_c == Constant.NULL) {
						neighbor[l_c][c] = Constant.NULL;
					}
				} else if (l_c == Constant.CAP0) {
					cap0_adj.get(c).remove((Integer) left);
					if (r_c != Constant.CAP0 && r_c != Constant.CAP1
							&& r_c != Constant.NULL) {
						cap0_adj.get(c).add(r_c);
						neighbor[r_c][c] = Constant.CAP0;
					} else if (r_c == Constant.CAP0) {
						cap0_adj.get(c).remove((Integer) right);
						path00[c]++;
					} else if (r_c == Constant.CAP1) {
						cap1_adj.get(c).remove((Integer) right);
						cycle_number++;
					} else if (r_c == Constant.NULL)
						;
				} else if (l_c == Constant.CAP1) {
					cap1_adj.get(c).remove((Integer) left);
					if (r_c != Constant.CAP0 && r_c != Constant.CAP1
							&& r_c != Constant.NULL) {
						cap1_adj.get(c).add(r_c);
						neighbor[r_c][c] = Constant.CAP1;
					} else if (r_c == Constant.CAP0) {
						cap0_adj.get(c).remove((Integer) right);
						cycle_number++;
					} else if (r_c == Constant.CAP1) {
						cap1_adj.get(c).remove((Integer) right);
						path11[c]++;
					} else if (r_c == Constant.NULL)
						;
				} else if (l_c == Constant.NULL) {
					if (r_c != Constant.CAP0 && r_c != Constant.CAP1
							&& r_c != Constant.NULL) {
						neighbor[r_c][c] = Constant.NULL;
					} else if (r_c == Constant.CAP0) {
						cap0_adj.get(c).remove((Integer) right);
					} else if (r_c == Constant.CAP1) {
						cap0_adj.get(c).remove((Integer) right);
					} else if (r_c == Constant.NULL)
						;
				}
				neighbor[left][c] = Constant.NULL;
				neighbor[right][c] = Constant.NULL;
				remains[map[left]] = remains[map[right]] = false;
			}
		} else if (left == Constant.CAP0 && right == Constant.CAP0) {
			for (int c = 0; c < 3; c++)
				path11[c]++;
			num_free_caps -= 2;
		} else if (left != Constant.CAP0 && right == Constant.CAP0) {
			adding_order[map[left]] = order;
			median_adj.reg_adj[map[left]] = Constant.CAP0;
			median_adj.cap_adj.add(map[left]);
			for (int c = 0; c < 3; c++) {
				int l_c = neighbor[left][c];
				if (l_c != Constant.CAP0 && l_c != Constant.CAP1
						&& l_c != Constant.NULL) {
					cap1_adj.get(c).add(l_c);
					neighbor[l_c][c] = Constant.CAP1;
				}
				if (l_c == Constant.CAP0) {
					cap0_adj.get(c).remove((Integer) left);
					cycle_number++;
				}
				if (l_c == Constant.CAP1) {
					cap1_adj.get(c).remove((Integer) left);
					path11[c]++;
				}
				if (l_c == Constant.NULL)
					;
				neighbor[left][c] = Constant.NULL;
				remains[map[left]] = false;
			}
			num_free_caps--;
		} else if (left == Constant.CAP0 && right != Constant.CAP0) {
			adding_order[map[right]] = order;
			median_adj.reg_adj[map[right]] = Constant.CAP0;
			median_adj.cap_adj.add(map[right]);
			for (int c = 0; c < 3; c++) {
				int r_c = neighbor[right][c];
				if (r_c != Constant.CAP0 && r_c != Constant.CAP1
						&& r_c != Constant.NULL) {
					cap1_adj.get(c).add(r_c);
					neighbor[r_c][c] = Constant.CAP1;
				}
				if (r_c == Constant.CAP0) {
					cap0_adj.get(c).remove((Integer) right);
					cycle_number++;
				}
				if (r_c == Constant.CAP1) {
					cap1_adj.get(c).remove((Integer) right);
					path11[c]++;
				}
				if (r_c == Constant.NULL)
					;
				neighbor[right][c] = Constant.NULL;
				remains[map[right]] = false;
			}
			num_free_caps--;
		}

		// check path00 and path11 to make cycle number right

		for (int c = 0; c < 3; c++) {
			if (path11[c] > 0 && path00[c] > 0) {
				int cyc = Math.min(path11[c], path00[c]);
				cycle_number += cyc;
				path11[c] -= cyc;
				path00[c] -= cyc;
			}
		}

		for (int i = 0; i < num_reg_vtx; i++)
			for (int c = 0; c < 3; c++) {
				if (neighbor[i][c] == i)
					System.out.printf("NULL should not exist! v=%d c=%d\n", i,
							c);
			}
	}

	public void shrink(ArrayList<Integer> black) {
		for (Integer b : black)
			if (b != Constant.CAP0 && b != Constant.CAP1)
				steps.add(map[b]);
			else
				steps.add(b);
		// for(Integer b:steps) System.out.print(b+" ");
		// System.out.println();
		int m = black.size() / 2; // the number of black edges, i.e. the size
		// of the subgraph

		// refactor the neighbor array and the cap array
		int cap_edges = 0;
		for (int i = 0; i < m; i++) {
			int l = black.get(2 * i); // the left vertex
			int r = black.get(2 * i + 1); // the right vertex
			if (l == r)
				continue;
			shrink(l, r);
		}

		relabel();
		if (confident)
			order++;

		// when we need null chromosomes
		if (num_free_caps < 0) {
			System.out.println("num_free_caps is smaller than 0");
			System.exit(1);
		}

		if (num_free_caps > 0 && num_reg_vtx == 0) {
			// System.out.println("number of chromosomes:"+(num_chr-num_free_caps/2));
			for (int color = 0; color < 3; color++) {
				if (num_free_caps / 2 != size) {
					System.out.println("size is note equal to num_free_caps/2");
					System.exit(1);
				}
				path11[color] += num_free_caps / 2;
				if (path00[color] > 0) {
					int min = Math.min(path00[color], path11[color]);
					if (path00[color] != path11[color]) {
						System.out.println("different number of paths");
						System.exit(1);
					}
					cycle_number += min;
					path00[color] -= min;
					path11[color] -= min;
				}
			}
			ArrayList<Integer> b = new ArrayList<Integer>(num_free_caps);
			for (int i = 0; i < num_free_caps; i++)
				b.add(Constant.CAP0);
			steps.addAll(b);
			size = num_free_caps = num_reg_vtx = 0;

		}

		get_bounds();
		int tmp[] = median_adj.chromosomeCount();
		num_plasmid = tmp[1];
		median_adj.num_chr = tmp[0];
	}

	public void relabel() {
		// we relabel two types of relationships
		// one for the mapping between original labels and new labels
		// one for the mapping between old labels and new labels
		// temporary arrays (old_new, new_old) are used

		int old_new[] = new int[num_reg_vtx];
		for (int i = 0; i < num_reg_vtx; i++)
			old_new[i] = i;

		int shift = 0;
		for (int i = 0; i < 2 * num_gene; i++) {
			if (!remains[i] && new_label[i] != Constant.NULL) {
				// new deleted vertex
				old_new[new_label[i]] = Constant.NULL;
				new_label[i] = Constant.NULL;
				shift++;
			} else if (remains[i]) {
				old_new[new_label[i]] -= shift;
				new_label[i] -= shift;
			}
		}
		num_reg_vtx -= shift;
		map = new int[num_reg_vtx];
		int new_old[] = new int[num_reg_vtx];
		for (int i = 0; i < 2 * num_gene; i++) {
			if (remains[i]) {
				map[new_label[i]] = i;
			}
		}
		for (int i = 0; i < num_reg_vtx + shift; i++)
			if (old_new[i] != Constant.NULL)
				new_old[old_new[i]] = i;

		int temp[][] = neighbor;

		neighbor = new int[num_reg_vtx][3];
		for (int i = 0; i < num_reg_vtx; i++) {
			int old_label = new_old[i];
			for (int color = 0; color < 3; color++) {
				int t = temp[old_label][color];

				if (t == Constant.NULL) {
					System.out.println("Strange!");
				}

				if (t == Constant.CAP0)
					neighbor[i][color] = Constant.CAP0;
				else if (t == Constant.CAP1)
					neighbor[i][color] = Constant.CAP1;
				else if (t == Constant.NULL)
					neighbor[i][color] = Constant.NULL;
				else
					neighbor[i][color] = old_new[t];
			}
		}

		// refactor the cap arrays
		ArrayList<ArrayList<Integer>> old_cap0_adj = cap0_adj;
		ArrayList<ArrayList<Integer>> old_cap1_adj = cap1_adj;
		cap0_adj = new ArrayList<ArrayList<Integer>>(3);
		cap1_adj = new ArrayList<ArrayList<Integer>>(3);
		for (int color = 0; color < 3; color++) {
			cap0_adj.add(new ArrayList<Integer>(old_cap0_adj.get(color).size()));
			cap1_adj.add(new ArrayList<Integer>(old_cap1_adj.get(color).size()));
		}

		for (int color = 0; color < 3; color++) {
			for (int s : old_cap0_adj.get(color)) {
				cap0_adj.get(color).add((Integer) old_new[s]);
			}
			for (int s : old_cap1_adj.get(color)) {
				cap1_adj.get(color).add((Integer) old_new[s]);
			}
		}

		for (int i = 0; i < num_reg_vtx; i++)
			for (int c = 0; c < 3; c++) {
				if (neighbor[i][c] == Constant.NULL)
					System.out.printf("NULL should not exist! v=%d c=%d\n", i,
							c);
			}
	}

	public int incident(int l, int r) {
		if (l >= Constant.CAP1 || r >= Constant.CAP1)
			return -1;
		for (int color = 0; color < 3; color++) {
			if (neighbor[l][color] == r)
				return color;
		}
		return -1;
	}

	public boolean is_connected(int l, int r) {
		if (l >= Constant.CAP1 || r >= Constant.CAP1)
			return false;
		for (int color = 0; color < 3; color++) {
			if (neighbor[l][color] == r)
				return true;
		}
		return false;
	}

	public int two_connected_by(int l, int r) {
		if (l == r)
			return -1;
		if (l >= Constant.CAP1 || r >= Constant.CAP1)
			return -1;
		for (int color = 0; color < 3; color++) {
			int lc = neighbor[l][color];
			if (lc >= Constant.CAP1)
				continue;
			for (int i = 1; i <= 2; i++) {
				int c = (color + i) % 3;
				if (lc == neighbor[l][c])
					return lc;
			}
		}
		return -1;
	}

	public void get_bounds() {
		int c01 = count_cycle(0, 1);
		int c02 = count_cycle(0, 2);
		int c12 = count_cycle(1, 2);
		upper_bound = cycle_number
				+ (int) Math.floor((3 * size + c01 + c02 + c12) / 2.0);
		int tc[] = new int[] { c01 + c02, c01 + c12, c02 + c12 };
		int tcm = tc[0], index = 0;
		for (int i = 1; i < 3; i++) {
			if (tc[i] > tcm) {
				tcm = tc[i];
				index = i;
			}
		}
		lower_bound = cycle_number + tcm;
	}

	public int count_cycle(int c1, int c2) {
		int cycles = 0;
		int pathaa = path00[c1], pathbb = path00[c2], pathab = 0, pathapap = path11[c1], pathbpbp = path11[c2], pathapbp = 0, pathaap = 0, pathbbp = 0, pathabp = 0, pathbap = 0;

		ArrayList<Integer> a = new ArrayList<Integer>(cap0_adj.get(c1));
		ArrayList<Integer> b = new ArrayList<Integer>(cap0_adj.get(c2));
		ArrayList<Integer> ap = new ArrayList<Integer>(cap1_adj.get(c1));
		ArrayList<Integer> bp = new ArrayList<Integer>(cap1_adj.get(c2));

		boolean unused[] = new boolean[num_reg_vtx];
		for (int i = 0; i < num_reg_vtx; i++)
			unused[i] = true;

		int left, right;
		while (!a.isEmpty()) {
			left = a.remove(a.size() - 1);
			unused[left] = false;
			while (true) {
				right = neighbor[left][c2];
				if (right == Constant.CAP0) {
					pathab++;
					b.remove((Integer) left);
					break;
				} else if (right == Constant.CAP1) {
					pathabp++;
					bp.remove((Integer) left);
					break;
				}
				unused[right] = false;

				left = neighbor[right][c1];
				if (left == Constant.CAP0) {
					pathaa++;
					a.remove((Integer) right);
					break;
				} else if (left == Constant.CAP1) {
					pathaap++;
					ap.remove((Integer) right);
					break;
				}
				unused[left] = false;
			}
		}

		while (!b.isEmpty()) {
			left = b.remove(b.size() - 1);
			unused[left] = false;
			while (true) {
				right = neighbor[left][c1];
				if (right == Constant.CAP1) {
					pathbap++;
					ap.remove((Integer) left);
					break;
				}
				unused[right] = false;

				left = neighbor[right][c2];
				if (left == Constant.CAP0) {
					pathbb++;
					b.remove((Integer) right);
					break;
				} else if (left == Constant.CAP1) {
					pathbbp++;
					bp.remove((Integer) right);
					break;
				}
				unused[left] = false;
			}
		}

		while (!ap.isEmpty()) {
			left = ap.remove(ap.size() - 1);
			unused[left] = false;
			while (true) {
				right = neighbor[left][c2];
				if (right == Constant.CAP1) {
					pathapbp++;
					bp.remove((Integer) left);
					break;
				}
				unused[right] = false;

				left = neighbor[right][c1];
				if (left == Constant.CAP1) {
					pathapap++;
					ap.remove((Integer) right);
					break;
				}
				unused[left] = false;
			}
		}

		while (!bp.isEmpty()) {
			left = bp.remove(bp.size() - 1);
			unused[left] = false;
			while (true) {
				right = neighbor[left][c1];
				unused[right] = false;

				left = neighbor[right][c2];
				if (left == Constant.CAP1) {
					pathbpbp++;
					bp.remove((Integer) right);
					break;
				}
				unused[left] = false;
			}
		}

		int start;
		for (int i = 0; i < num_reg_vtx; i++) {
			if (!unused[i])
				continue;
			start = left = i;
			do {
				right = neighbor[left][c1];

				if (left == Constant.CAP0)
					System.out.println("left is a cap");
				else if (left == Constant.CAP1)
					System.out.println("left is a cap1");
				else if (left == Constant.NULL)
					System.out.println("left is a null");
				if (right == Constant.CAP0)
					System.out.println("right is a cap");
				else if (right == Constant.CAP1)
					System.out.println("right is a cap1");
				else if (right == Constant.NULL)
					System.out.println("right is a null");

				unused[left] = unused[right] = false;
				left = neighbor[right][c2];
			} while (left != start);
			cycles++;
		}

		// System.out.println("in calculating cycles");
		// System.out.println("pathaa:"+pathaa+"\tpathbb:"+pathbb+"\tpathab:"+pathab);
		// System.out.println("pathapap:"+pathapap+"\tpathbpbp:"+pathbpbp+"\tpathapbp:"+pathapbp);
		// System.out.println("pathaap:"+pathaap+"\tpathbbp:"+pathbbp+"\tpathabp:"+pathabp+"\tpathbap:"+pathbap);

		// if (pathab >= num_free_caps) {
		// cycles += num_free_caps;
		// pathab -= num_free_caps;
		// } else {
		// cycles += pathab + (num_free_caps - pathab) / 2;
		// pathaa -= (num_free_caps - pathab) / 2;
		// pathbb -= (num_free_caps - pathab) / 2;
		// pathab = 0;
		// }

		pathapbp += num_free_caps;

		// if((pathabp + pathbap+ pathab + pathapbp)%2!=0) {
		// System.out.println("error in calculating cycle numbers"+"\tnum_free_caps:"+num_free_caps);
		// System.exit(1);
		// }

		return cycles + pathaap + pathbbp
				+ (pathabp + pathbap + pathab + pathapbp) / 2
				+ (pathaa < pathapap ? pathaa : pathapap)
				+ (pathbb < pathbpbp ? pathbb : pathbpbp);
	}

	public int getCycle_number() {
		return cycle_number;
	}

}
