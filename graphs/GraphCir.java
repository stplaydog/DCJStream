package graphs;

import order.GeneOrder;
import tools.Params;
import gasts.*;

public class GraphCir extends Graph {

	public void init(GeneOrder order) {
		// TODO Auto-generated constructor stub
		this.size = order.gene_num * 2;
		this.gene_num = order.gene_num;
		this.node_num = gene_num * 2;
		this.edge_num = gene_num;
		adj_mat = new int[node_num][3];
		vet_rank = new int[node_num][3];
		vet_cynum = new int[node_num][3];
		this.cynum = new int[node_num][3];
		this.cyrank = new int[node_num][3];
		check = new boolean[node_num];
		unused = new boolean[node_num];
		cycle = new int[3];
		cycle_number = 0;
		c = new int[3];
		lower_bound = 0;
		upper_bound = 0;
		footprint = new int[this.node_num];
		major_tmp = new int[this.node_num];
		this.idx_ft = 0;
		this.idx_tmp = 0;
		this.gene_order_to_graph(order);
		this.get_bounds();
		// initialize the adj
		ft_before_idx = 0;
		map = new int[size];
		ft_before_rename = new int[size];
		median_adj = new Adjacency(gene_num, 1, "");
		median_adj.reg_adj[0] = median_adj.reg_adj[gene_num * 2 - 1];
		median_adj.reg_adj[gene_num * 2 - 1] = median_adj.reg_adj[0];
		for (int i = 0; i < gene_num - 1; i++) {
			median_adj.reg_adj[i * 2 + 1] = median_adj.reg_adj[(i + 1) * 2];
			median_adj.reg_adj[(i + 1) * 2] = median_adj.reg_adj[i * 2 + 1];
		}
	}

	public void init(Adjacency adj1, Adjacency adj2, Adjacency adj3) {
		// TODO Auto-generated constructor stub
		this.size = adj1.num_gene * 2;
		this.gene_num = adj1.num_gene;
		this.node_num = gene_num * 2;
		this.edge_num = gene_num;
		adj_mat = new int[node_num][3];
		vet_rank = new int[node_num][3];
		vet_cynum = new int[node_num][3];
		this.cynum = new int[node_num][3];
		this.cyrank = new int[node_num][3];
		check = new boolean[node_num];
		unused = new boolean[node_num];
		cycle = new int[3];
		cycle_number = 0;
		c = new int[3];
		lower_bound = 0;
		upper_bound = 0;
		footprint = new int[this.node_num];
		major_tmp = new int[this.node_num];
		this.idx_ft = 0;
		this.idx_tmp = 0;
		this.gene_order_to_graph(adj1, adj2, adj3);
		this.get_bounds();
		// initialize the adj
		ft_before_idx = 0;
		map = new int[size];
		ft_before_rename = new int[size];
		median_adj = new Adjacency(gene_num, 1, "");
		median_adj.reg_adj[0] = median_adj.reg_adj[gene_num * 2 - 1];
		median_adj.reg_adj[gene_num * 2 - 1] = median_adj.reg_adj[0];
		for (int i = 0; i < gene_num - 1; i++) {
			median_adj.reg_adj[i * 2 + 1] = median_adj.reg_adj[(i + 1) * 2];
			median_adj.reg_adj[(i + 1) * 2] = median_adj.reg_adj[i * 2 + 1];
		}
	}

	@Override
	public void init(Graph g) {
		// TODO Auto-generated constructor stub
		this.size = g.gene_num * 2;
		this.gene_num = g.gene_num;
		this.node_num = gene_num * 2;
		this.edge_num = g.edge_num;
		adj_mat = new int[node_num][3];
		for (int i = 0; i < 3; i++)
			for (int j = 0; j < node_num; j++)
				this.adj_mat[j][i] = g.adj_mat[j][i];
		vet_rank = new int[node_num][3];
		vet_cynum = new int[node_num][3];
		this.cynum = new int[node_num][3];
		this.cyrank = new int[node_num][3];
		unused = new boolean[node_num];
		check = new boolean[node_num];
		for (int j = 0; j < node_num; j++)
			this.check[j] = g.check[j];
		cycle = new int[3];
		cycle_number = g.cycle_number;
		c = new int[3];
		lower_bound = g.lower_bound;
		upper_bound = g.upper_bound;
		this.idx_ft = g.idx_ft;
		footprint = new int[this.node_num];
		for (int i = 0; i < g.idx_ft; i++)
			this.footprint[i] = g.footprint[i];
		major_tmp = new int[this.node_num * 2];
		this.idx_tmp = g.idx_tmp;
		this.v_num = g.v_num;
	}

	public void gene_order_to_graph(Adjacency adj1, Adjacency adj2,
			Adjacency adj3) {
		for (int i = 0; i < node_num; i++) {
			adj_mat[i][0] = adj1.reg_adj[i];
			adj_mat[i][1] = adj2.reg_adj[i];
			adj_mat[i][2] = adj3.reg_adj[i];
		}

		cycle_number = 0;
		for (int i = 0; i < node_num; i++)
			check[i] = true;

	}

	public void gene_order_to_graph(GeneOrder o) {
		int color, gene;
		int i;
		int order, order_forward, order_backward;
		int forward_pos, backward_pos;
		int forward_pt, backward_pt;
		// double direction edge

		for (color = 0; color < 3; color++) {
			for (gene = 0; gene < gene_num; gene++) {
				// (forward_pt)<-(forward_pos) (backward_pos)->(backward_pt)
				order = o.gene_order[gene][color];
				// forward
				forward_pos = order > 0 ? (Math.abs(order) - 1) * 2 : (Math
						.abs(order) - 1) * 2 + 1;
				order_forward = o.gene_order[(gene - 1) >= 0 ? (gene - 1)
						: (gene_num - 1)][color];
				forward_pt = order_forward > 0 ? (Math.abs(order_forward) - 1) * 2 + 1
						: (Math.abs(order_forward) - 1) * 2;
				adj_mat[forward_pos][color] = forward_pt;
				// backward
				backward_pos = order > 0 ? (Math.abs(order) - 1) * 2 + 1
						: (Math.abs(order) - 1) * 2;
				order_backward = o.gene_order[(gene + 1) < gene_num ? (gene + 1)
						: (0)][color];
				backward_pt = order_backward > 0 ? (Math.abs(order_backward) - 1) * 2
						: (Math.abs(order_backward) - 1) * 2 + 1;
				adj_mat[backward_pos][color] = backward_pt;
			}
		}
		cycle_number = 0;
		for (i = 0; i < node_num; i++)
			check[i] = true;
	}

	public void shrink(int[] major, int start, int end) {
		int i, j;
		int left, right;
		boolean valid[] = new boolean[node_num];

		for (i = 0; i < node_num; i++) {
			if (check[i] == true)
				valid[i] = true;
			else
				valid[i] = false;
		}

		for (i = start; i < end; i++) {
			valid[major[i]] = false;
			check[major[i]] = false;
		}

		// the real shrinking part
		for (i = start; i < end; i += 2) {
			left = major[i];
			right = major[i + 1];
			// color1
			if (adj_mat[left][0] == right) {
				cycle_number++;
			} else {
				adj_mat[adj_mat[left][0]][0] = adj_mat[right][0];
				adj_mat[adj_mat[right][0]][0] = adj_mat[left][0];
			}
			// color2
			if (adj_mat[left][1] == right) {
				cycle_number++;
			} else {
				adj_mat[adj_mat[left][1]][1] = adj_mat[right][1];
				adj_mat[adj_mat[right][1]][1] = adj_mat[left][1];
			}
			// color3
			if (adj_mat[left][2] == right) {
				cycle_number++;
			} else {
				adj_mat[adj_mat[left][2]][2] = adj_mat[right][2];
				adj_mat[adj_mat[right][2]][2] = adj_mat[left][2];
			}
			check[left] = false;
			check[right] = false;
		}

		for (j = start; j < end; j++)
			this.add_ft(major[j]);

		edge_num -= (end - start) / 2;
	}

	public void expand(int[] major, int start, int end) {
		int i;
		int left, right;

		for (i = start; i < end; i++) {
			check[major[i]] = true;
		}

		// the real shrinking part
		for (i = end - 1; i >= start; i -= 2) {
			left = major[i - 1];
			right = major[i];
			// color1
			if (adj_mat[left][0] == right) {
				cycle_number--;
			} else {
				adj_mat[adj_mat[left][0]][0] = left;
				adj_mat[adj_mat[right][0]][0] = right;
			}
			// color2
			if (adj_mat[left][1] == right) {
				cycle_number--;
			} else {
				adj_mat[adj_mat[left][1]][1] = left;
				adj_mat[adj_mat[right][1]][1] = right;
			}
			// color3
			if (adj_mat[left][2] == right) {
				cycle_number--;
			} else {
				adj_mat[adj_mat[left][2]][2] = left;
				adj_mat[adj_mat[right][2]][2] = right;
			}
			check[left] = true;
			check[right] = true;
		}
		for (i = start; i < end; i++)
			this.mv_back_ft();
		edge_num += (end - start) / 2;
	}

	public void rename() {
		int new_num_vet = 0;
		int idx = 0;
		map = new int[node_num];
		for (int i = 0; i < gene_num * 2; i++) {
			if (check[i] == true) {
				map[i] = idx;
				idx++;
				new_num_vet++;
			} else
				map[i] = -1;
		}
		idx = 0;
		for (int i = 0; i < gene_num * 2; i++) {
			if (check[i] == true) {
				int c1 = adj_mat[i][0];
				int c2 = adj_mat[i][1];
				int c3 = adj_mat[i][2];
				adj_mat[idx][0] = map[c1];
				adj_mat[idx][1] = map[c2];
				adj_mat[idx][2] = map[c3];
				idx++;
			}
		}

		for (int i = 0; i < new_num_vet; i++)
			check[i] = true;

		node_num = new_num_vet;
		gene_num = new_num_vet / 2;
		edge_num = new_num_vet / 2;
		v_num = new_num_vet;
	}

	@Override
	public void init(AdjNode n1, AdjNode n2, AdjNode n3) {
		// TODO Auto-generated method stub

		// TODO Auto-generated constructor stub
		this.size = n1.adj.num_gene * 2;
		this.gene_num = n1.adj.num_gene;
		this.node_num = gene_num * 2;
		this.edge_num = n1.adj.num_gene;
		adj_mat = new int[node_num][3];
		for (int j = 0; j < node_num; j++) {
			this.adj_mat[j][0] = n1.adj.reg_adj[j];
			this.adj_mat[j][1] = n2.adj.reg_adj[j];
			this.adj_mat[j][2] = n3.adj.reg_adj[j];
		}
		vet_rank = new int[node_num][3];
		vet_cynum = new int[node_num][3];
		this.cynum = new int[node_num][3];
		this.cyrank = new int[node_num][3];
		unused = new boolean[node_num];
		check = new boolean[node_num];
		for (int j = 0; j < node_num; j++)
			this.check[j] = true;
		cycle = new int[3];
		cycle_number = 0;
		c = new int[3];
		lower_bound = 0;
		upper_bound = 0;
		this.idx_ft = 0;
		footprint = new int[this.node_num];
		for (int i = 0; i < this.node_num; i++)
			this.footprint[i] = 0;
		major_tmp = new int[this.node_num * 2];
		this.idx_tmp = 0;
		this.v_num = this.node_num;

	}

	@Override
	public Adjacency toMedianAdj() {
		// copy the renamed vertices
		for (int i = 0; i < ft_before_idx; i += 2) {
			median_adj.reg_adj[ft_before_rename[i]] = ft_before_rename[i + 1];
			median_adj.reg_adj[ft_before_rename[i + 1]] = ft_before_rename[i];
		}
		// copy the current ft vertices
		for (int i = 0; i < this.idx_ft; i += 2) {
			median_adj.reg_adj[footprint[i]] = footprint[i + 1];
			median_adj.reg_adj[footprint[i + 1]] = footprint[i];
		}
		// greedy add nodes
		int from_adj = -1;
		int to_adj = -1;
		while ((from_adj = this.next_median_adj()) != -1) {
			to_adj = this.next_median_adj();
			median_adj.reg_adj[from_adj]=to_adj;
			median_adj.reg_adj[to_adj]=from_adj;
		}
		return this.median_adj;
	}

	public void get_bounds() {
		int lower_ind = -1;
		lower_ind = -1;
		count_cycle(0, 1);
		count_cycle(0, 2);
		count_cycle(1, 2);

		if (c[0] <= c[1] && c[0] <= c[2])
			lower_ind = 0;
		else if (c[1] <= c[2] && c[1] <= c[2])
			lower_ind = 1;
		else if (c[2] <= c[0] && c[2] <= c[1])
			lower_ind = 2;
		// up = cycle_number + (3 * edge_num + c[0] + c[1] + c[2]) / 2;
		// dis_to_cor();
		// upper_bound = cycle_number +3 * edge_num - cor_to_lb();
		upper_bound = cycle_number
				+ (int) Math.floor((3 * edge_num + c[0] + c[1] + c[2]) / 2);
		lower_bound = cycle_number + edge_num + c[0] + c[1] + c[2]
				- c[lower_ind];
	}

	public void count_cycle(int c1, int c2) {
		int i, cycles;
		int start, left, right;
		int co;
		boolean unused[] = new boolean[node_num];

		if (c1 == 0 && c2 == 1)
			co = 2;
		else if (c1 == 0 && c2 == 2)
			co = 1;
		else
			co = 0;

		cycles = 0;
		for (i = 0; i < node_num; i++) {
			if (check[i] == true)
				unused[i] = true;
			else
				unused[i] = false;
		}

		for (i = 0; i < node_num; i++) {
			if (unused[i] != true)
				continue;
			start = left = i;
			do {
				right = adj_mat[left][c1];
				unused[left] = unused[right] = false;
				left = adj_mat[right][c2];
			} while (left != start);
			cycles++;
		}
		c[co] = cycles;
	}

	public void get_rank_cynum(int c1, int c2) {
		int i, cycles;
		int start, left, right;

		int co;
		int rank;

		cycles = 0;

		boolean unused[] = new boolean[node_num];

		if (c1 == 0 && c2 == 1)
			co = 2;
		else if (c1 == 0 && c2 == 2)
			co = 1;
		else
			co = 0;

		for (i = 0; i < node_num; i++) {
			if (check[i] == true)
				unused[i] = true;
			else
				unused[i] = false;
		}

		for (i = 0; i < node_num; i++) {
			if (unused[i] != true)
				continue;
			start = left = i;
			rank = 0;
			do {
				right = adj_mat[left][c1];
				unused[left] = unused[right] = false;
				vet_cynum[left][co] = cycles;
				vet_rank[left][co] = rank;
				vet_cynum[right][co] = cycles;
				vet_rank[right][co] = rank + 1;
				rank += 2;
				left = adj_mat[right][c2];
			} while (left != start);
			cycles++;
		}
		c[co] = cycles;
	}

	public void get_bounds_linear(int v1, int v2) {
		int lower_ind = -1;
		lower_ind = -1;
		count_cycle_linear(v1, v2, 0, 1);
		count_cycle_linear(v1, v2, 0, 2);
		count_cycle_linear(v1, v2, 1, 2);

		if (c[0] <= c[1] && c[0] <= c[2])
			lower_ind = 0;
		else if (c[1] <= c[2] && c[1] <= c[2])
			lower_ind = 1;
		else if (c[2] <= c[0] && c[2] <= c[1])
			lower_ind = 2;
		// up = cycle_number + (3 * edge_num + c[0] + c[1] + c[2]) / 2;
		// dis_to_cor();
		// upper_bound = cycle_number + 3 * edge_num-cor_to_lb();
		// upper_bound = cycle_number + (3 * edge_num + c[0] + c[1] + c[2]) / 2;
		upper_bound = cycle_number
				+ (int) Math.floor((3 * edge_num + c[0] + c[1] + c[2]) / 2);
		lower_bound = cycle_number + edge_num + c[0] + c[1] + c[2]
				- c[lower_ind];
	}

	public void count_cycle_linear(int v1, int v2, int c1, int c2) {
		int color;
		int x, y, w, z;
		int r_x, r_y, r_w, r_z;

		if (c1 == 0 && c2 == 1)
			color = 2;
		else if (c1 == 0 && c2 == 2)
			color = 1;
		else
			color = 0;

		if (vet_cynum[v1][color] != vet_cynum[v2][color])
			c[color]--;
		else {
			x = adj_mat[v1][c1];
			y = adj_mat[v1][c2];
			w = adj_mat[v2][c1];
			z = adj_mat[v2][c2];
			r_x = vet_rank[x][color];
			r_y = vet_rank[y][color];
			r_w = vet_rank[w][color];
			r_z = vet_rank[z][color];

			if (r_x == r_z || r_y == r_w)
				return;
			else if (x == v2 || y == v2 || w == v1 || z == v1)
				return;

			int left, right;
			int start = left = x;
			do {
				right = adj_mat[left][c2];
				left = adj_mat[right][c1];
				if (left == y || left == z)
					return;
			} while (left != start);
			c[color]++;

		}
	}

	public void get_start(int c1, int c2) {
		int i, cycles;
		int start, left, right;
		int color;

		if (c1 == 0 && c2 == 1)
			color = 2;
		else if (c1 == 0 && c2 == 2)
			color = 1;
		else
			color = 0;

		int rank;

		cycles = 0;

		for (i = 0; i < node_num; i++) {
			unused[i] = check[i];
		}

		for (i = 0; i < node_num; i++) {
			if (unused[i] != true)
				continue;
			start = left = i;
			rank = 0;
			do {
				right = adj_mat[left][c1];
				unused[left] = unused[right] = false;
				cynum[left][color] = cycles;
				cynum[right][color] = cycles;
				rank += 2;
				left = adj_mat[right][c2];
			} while (left != start);
			cyrank[cycles][color] = rank;
			cycles++;
		}
	}

	public String toString() {
		StringBuffer buffer = new StringBuffer();
		int new_num_vet = 0;
		int idx = 0;
		int map[] = new int[node_num];
		for (int i = 0; i < gene_num * 2; i++) {
			if (check[i] == true) {
				map[i] = idx;
				idx++;
				new_num_vet++;
			} else
				map[i] = -1;
		}
		idx = 0;
		int new_adj[][] = new int[new_num_vet][3];
		for (int i = 0; i < gene_num * 2; i++) {
			if (check[i] == true) {
				int c1 = adj_mat[i][0];
				int c2 = adj_mat[i][1];
				int c3 = adj_mat[i][2];
				new_adj[idx][0] = map[c1];
				new_adj[idx][1] = map[c2];
				new_adj[idx][2] = map[c3];
				buffer.append("[" + i + ":" + new_adj[idx][0] + ","
						+ new_adj[idx][1] + "," + new_adj[idx][2] + "]" + ",");
				buffer.append("[" + i + ":" + c1 + "," + c2 + "," + c3 + "]"
						+ "\n");
				idx++;
			}
		}

		buffer.append("\nub:" + upper_bound + " lb:" + lower_bound + "\n");
		buffer.append(this.footprint.toString() + "\n");
		return buffer.toString();
	}

	public String toHash() {
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < this.idx_ft; i++)
			buffer.append(this.footprint[i]);
		return buffer.toString();
	}

	public void toOrder(Params p, GeneOrder order) {
		for (int color = 0; color < 3; color++) {
			for (int i = 0; i < node_num; i++)
				this.unused[i] = true;

			int idx = 0;
			int gene;
			int current = -1, next;
			for (short i = 0; i < node_num; i++) {
				if (!unused[i])
					continue;
				current = i;
				do {
					unused[current] = false;
					if (current % 2 == 0) {
						gene = current / 2 + 1;
						next = (short) (current + 1);
					} else {
						gene = (-1) * (current / 2 + 1);
						next = (short) (current - 1);
					}
					order.gene_order[idx][color] = gene;
					idx++;
					unused[next] = false;
					current = adj_mat[next][color];
				} while (current != i);
			}
		}
		this.gene_order_to_graph(order);
	}
}
