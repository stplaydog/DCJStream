package graphs;

import gasts.AdjNode;
import gasts.Adjacency;

import java.util.ArrayList;

import order.*;
import tools.*;

public abstract class Graph {
	public boolean is_heu = false;
	public int node_num; // how many nodes
	public int gene_num; // how many genes
	public int edge_num; // how many edges
	public int chr_num; // how many chromosomes
	public int v_num;
	public int adj_mat[][]; // adjacency matrix, example: +1->0 -1->1 +2->2
							// -2->3
	public int vet_rank[][];
	public int vet_cynum[][];

	public CList cap_vet[];
	public CList cup_vet[];

	int cynum[][];
	int cyrank[][];

	public boolean check[]; // check which vertex is available

	public int cycle[]; // cycle number by
	public int cycle_number; // how many total cycles
	public int c[];

	public int upper_bound; // upper bound
	public int lower_bound; // lower bound

	public int footprint[];
	public int idx_ft;

	public int major_tmp[];
	public int idx_tmp;

	CList a;
	CList b;
	CList ap;
	CList bp;

	int path00[];
	int path11[];
	int cycle00[];
	int path_final[];
	int cycle_ambi[];

	int size;
	int chr_num_final;

	boolean unused[];

	int pa_count[][];
	int vet_type[][];
	int vet_panum[][];
	int spec_path[][];
	int pa_path[][];

	int ft_before_idx;
	public int map[];
	public int map_o[];
	public int ft_before_rename[];
	public Adjacency median_adj;

	public abstract void init(GeneOrder order);

	public abstract void init(AdjNode n1, AdjNode n2, AdjNode n3);

	public abstract void init(Graph g);

	public abstract void gene_order_to_graph(GeneOrder o);

	public abstract void shrink(int[] major, int start, int end);

	public abstract void expand(int[] major, int start, int end);

	public abstract void rename();

	public abstract void get_bounds();

	public abstract void count_cycle(int c1, int c2);

	public abstract void get_rank_cynum(int c1, int c2);

	public abstract void get_bounds_linear(int v1, int v2);

	public abstract void count_cycle_linear(int v1, int v2, int c1, int c2);

	public abstract void get_start(int c1, int c2);

	public abstract String toHash();

	public abstract void toOrder(Params p, GeneOrder order);

	public abstract Adjacency toMedianAdj();

	public boolean check_mute_violate(int v1, int v2) {
		if (v1 > v2) {
			if (v1 - v2 == 1 && v1 % 2 == 1)
				return true;
		} else {
			if (v2 - v1 == 1 && v2 % 2 == 1)
				return true;
		}
		return false;
	}

	public void mute(double mute_rate, Params p) {
		// TODO Auto-generated method stub
		int cut1;
		int cut2;
		int point1;
		int point2;
		int is_cir1 = 0;
		int is_cir2 = 0;
		int tmp_times = (int) (mute_rate * (double) this.gene_num);
		int num_times = 0;
		for (int color = 0; color < 3; color++) {
			num_times = tmp_times;
			while (num_times-- > 0) {
				boolean changed = false;
				// System.out.println(num_times);
				int min = 0;
				int max = this.size - 1;
				cut1 = min + (int) (Math.random() * ((max - min) + 1));
				while (true) {
					cut2 = min + (int) (Math.random() * ((max - min) + 1));
					if (cut2 == cut1 || cut1 == adj_mat[cut2][color])
						continue;
					else if (adj_mat[cut1][color] == Const.CAP
							&& adj_mat[cut2][color] == Const.CAP)
						continue;
					else {
						point1 = adj_mat[cut1][color];
						point2 = adj_mat[cut2][color];
						if (p.type.equals("lin")) {
							is_cir1 = is_generate_cycle(cut1, cut2, color);
							is_cir2 = is_generate_cycle(cut2, cut1, color);
						}
						break;
					}
				}
				if (check_mute_violate(cut1, cut2)
						|| check_mute_violate(point1, point2))
					continue;
				// two different edge and they both not point to cap
				if (adj_mat[cut1][color] != Const.CAP
						&& adj_mat[cut2][color] != Const.CAP) {
					if (Math.random() < 0.5
							&& (is_cir1 == 0 || is_cir1 != is_cir2)) {
						adj_mat[cut1][color] = cut2;
						adj_mat[cut2][color] = cut1;
						adj_mat[point1][color] = point2;
						adj_mat[point2][color] = point1;
						changed = true;
					} else if ((is_cir1 == 0 || is_cir1 == is_cir2)) {
						adj_mat[cut1][color] = point2;
						adj_mat[cut2][color] = point1;
						adj_mat[point1][color] = cut2;
						adj_mat[point2][color] = cut1;
						changed = true;
					}

				} else if (adj_mat[cut1][color] == Const.CAP
						&& adj_mat[cut2][color] != Const.CAP) {
					if (Math.random() < 0.5
							&& (is_cir1 == 0 || is_cir1 != is_cir2)) {
						adj_mat[cut1][color] = cut2;
						adj_mat[cut2][color] = cut1;
						cap_vet[color].remove(cut1);
						cap_vet[color].add(point2);
						adj_mat[point2][color] = Const.CAP;
						changed = true;
					} else if (is_cir1 == 0 || is_cir1 == is_cir2) {
						adj_mat[cut1][color] = point2;
						adj_mat[point2][color] = cut1;
						cap_vet[color].remove(cut1);
						cap_vet[color].add(cut2);
						adj_mat[cut2][color] = Const.CAP;
						changed = true;
					}
				} else if (adj_mat[cut1][color] != Const.CAP
						&& point2 == Const.CAP) {
					if (Math.random() < 0.5
							&& (is_cir1 == 0 || is_cir1 != is_cir2)) {
						adj_mat[cut1][color] = cut2;
						adj_mat[cut2][color] = cut1;
						cap_vet[color].remove(cut2);
						cap_vet[color].add(point1);
						adj_mat[point1][color] = Const.CAP;
						changed = true;
					} else if (is_cir1 == 0 || is_cir1 == is_cir2) {
						adj_mat[cut2][color] = point1;
						adj_mat[point1][color] = cut2;
						cap_vet[color].remove(cut2);
						cap_vet[color].add(cut1);
						adj_mat[cut1][color] = Const.CAP;
						changed = true;
					}
					if (changed == false)
						num_times++;
				}
			}

		}
		this.get_bounds();
	}

	public int connected_by(int l, int r) {
		int color;
		for (color = 0; color < 3; color++) {
			if (adj_mat[l][color] == r)
				return color;
		}
		return -1;
	}

	public int two_connected_by(int l, int r) {
		int i, color;
		int lc, c;
		if (l == r)
			return -1;
		for (color = 0; color < 3; color++) {
			lc = adj_mat[l][color];
			for (i = 1; i <= 2; i++) {
				c = (color + i) % 3;
				if (lc == adj_mat[r][c])
					return lc;
			}
		}
		return -1;
	}

	public boolean is_connected(int l, int r) {
		if (connected_by(l, r) == -1)
			return false;
		else
			return true;
	}

	public int incident(int l, int r) {
		if (l >= Const.CUP || r >= Const.CUP)
			return -1;
		for (int color = 0; color < 3; color++) {
			if (adj_mat[l][color] == r)
				return color;
		}
		return -1;
	}

	public void graph_vis_two(int l, int r) {
		int disable = -1;
		if (l == 0 && r == 1)
			disable = 2;
		else if (l == 0 && r == 2)
			disable = 1;
		else if (l == 1 && r == 2)
			disable = 0;

		System.out.print("graph G {\n");
		for (int i = 0; i < gene_num * 2; i++) {
			if (check[i] == true) {
				if (adj_mat[i][0] > i && disable != 0)
					System.out.printf("%d -- %d [color=red];\n", i,
							adj_mat[i][0]);
				if (adj_mat[i][1] > i && disable != 1)
					System.out.printf("%d -- %d [color=blue];\n", i,
							adj_mat[i][1]);
				if (adj_mat[i][2] > i && disable != 2)
					System.out.printf("%d -- %d [color=green];\n", i,
							adj_mat[i][2]);
			}
		}
		System.out.print("}\n");
	}

	public void graph_vis_one(int c) {
		int enable = c;

		System.out.print("graph G {\n");
		for (int i = 0; i < gene_num * 2; i++) {
			if (i % 2 == 0)
				System.out.printf("%d -- %d [color=black];\n", i, i + 1);
			if (check[i] == true) {
				if (adj_mat[i][0] > i && enable == 0)
					System.out.printf("%d -- %d [color=red];\n", i,
							adj_mat[i][0]);
				if (adj_mat[i][1] > i && enable == 1)
					System.out.printf("%d -- %d [color=blue];\n", i,
							adj_mat[i][1]);
				if (adj_mat[i][2] > i && enable == 2)
					System.out.printf("%d -- %d [color=green];\n", i,
							adj_mat[i][2]);
			}
		}
		System.out.print("}\n");
	}

	public void graph_vis_median() {
		System.out.print("graph G {\n");
		for (int i = 0; i < this.median_adj.num_gene * 2; i++) {
			if (i < this.median_adj.reg_adj[i])
				System.out.printf("%d -- %d [color=black];\n", i,
						this.median_adj.reg_adj[i]);
		}
		System.out.print("}\n");
	}

	public int next_median_adj() {
		for (int i = 0; i < this.median_adj.num_gene * 2; i++)
			if (this.median_adj.reg_adj[i] == Integer.MAX_VALUE) {
				this.median_adj.reg_adj[i] = -1;
				return i;
			}
		return -1;
	}

	public void graph_vis_all() {
		System.out.print("graph G {\n");
		for (int i = 0; i < gene_num * 2; i++) {
			if (check[i] == true) {
				if (adj_mat[i][0] > i)
					System.out.printf("%d -- %d [color=red];\n", i,
							adj_mat[i][0]);
				if (adj_mat[i][1] > i)
					System.out.printf("%d -- %d [color=blue];\n", i,
							adj_mat[i][1]);
				if (adj_mat[i][2] > i)
					System.out.printf("%d -- %d [color=green];\n", i,
							adj_mat[i][2]);
			}
		}
		System.out.print("}\n");
	}

	public void clean_ft() {
		// copy to a temp location
		for (int i = 0; i < this.idx_ft; i++, ft_before_idx++) {
			this.ft_before_rename[ft_before_idx] = this.footprint[i];
		}
		this.idx_ft = 0;
	}

	public void clean_major_tmp() {
		this.idx_tmp = 0;
	}

	public void add_ft(int v) {
		this.footprint[this.idx_ft] = v;
		this.idx_ft++;
	}

	public void mv_back_ft() {
		this.idx_ft--;
	}

	public int getIdx_ft() {
		return idx_ft;
	}

	public void setIdx_ft(int idx_ft) {
		this.idx_ft = idx_ft;
	}

	public int getIdx_tmp() {
		return idx_tmp;
	}

	public void setIdx_tmp(int idx_tmp) {
		this.idx_tmp = idx_tmp;
	}

	public int is_generate_cycle(int cut1, int cut2, int color) {
		int current;
		int next = cut1;
		while (next != Const.CAP) {
			current = next;
			if (current % 2 == 0) {
				next = (short) (current + 1);
			} else {
				next = (short) (current - 1);
			}
			if (next == cut2)
				return 1;
			current = next;
			next = adj_mat[current][color];
			if (next == cut2)
				return 2;
		}

		next = adj_mat[cut1][color];
		if (next == cut2)
			return 2;
		while (next != Const.CAP) {
			current = next;
			if (current % 2 == 0) {
				next = (short) (current + 1);
			} else {
				next = (short) (current - 1);
			}
			if (next == cut2)
				return 1;
			current = next;
			next = adj_mat[current][color];
			if (next == cut2)
				return 2;
		}
		return 0;
	}

	public int cal_up_start() {
		int sum = 0;
		int min_sum = 100000;
		int result = 0;
		for (int i = 0; i < this.v_num; i++) {
			if (check[i]) {
				sum = cyrank[cynum[i][0]][0] + cyrank[cynum[i][1]][1]
						+ cyrank[cynum[i][2]][2];
				if (sum < min_sum) {
					result = i;
					min_sum = sum;
				}
			}
		}
		return result;
	}

	public int cal_low_start() {
		int sum = 0;
		int max_sum = 0;
		int result = 0;
		for (int i = 0; i < this.v_num; i++) {
			if (check[i]) {
				sum = cyrank[cynum[i][0]][0] + cyrank[cynum[i][1]][1]
						+ cyrank[cynum[i][2]][2];
				if (sum > max_sum) {
					result = i;
					max_sum = sum;
				}
			}
		}
		return result;
	}

	public boolean check_f()
	{
		if(this.idx_ft!=8)
			return false;
		if(this.footprint[0]==8 &&
				this.footprint[1]==69 &&
				this.footprint[2]==55 &&
				this.footprint[3]==61 &&
				this.footprint[4]==3 &&
				this.footprint[5]==31 &&
				this.footprint[6]==53 &&
				this.footprint[7]==25 )
			return true;
		return false;
	}
}
