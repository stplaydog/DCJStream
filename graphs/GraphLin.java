package graphs;

import gasts.AdjNode;
import gasts.Adjacency;
import gasts.Constant;

import java.util.ArrayList;

import order.GeneOrder;
import tools.Const;
import tools.Params;

public class GraphLin extends Graph {

	@Override
	public void init(GeneOrder order) {
		// TODO Auto-generated constructor stub
		this.gene_num = order.gene_num;
		this.v_num = order.gene_num * 2;
		this.node_num = v_num;
		this.edge_num = gene_num;
		this.chr_num = order.chr_num[0] * 2;
		this.path00 = new int[3];
		this.path11 = new int[3];
		this.cycle00 = new int[3];
		this.cycle_ambi = new int[3];
		this.path_final = new int[3];
		adj_mat = new int[node_num][3];
		vet_rank = new int[node_num][3];
		vet_cynum = new int[node_num][3];
		check = new boolean[node_num];

		pa_count = new int[node_num][3];
		vet_type = new int[node_num][3];
		vet_panum = new int[node_num][3];
		spec_path = new int[node_num][3];
		pa_path = new int[node_num][3];

		this.cynum = new int[node_num][3];
		this.cyrank = new int[node_num][3];

		cycle = new int[3];
		cycle_number = 0;
		c = new int[3];
		lower_bound = 0;
		upper_bound = 0;
		footprint = new int[this.gene_num * 4];
		major_tmp = new int[this.node_num];
		this.unused = new boolean[this.node_num];
		this.idx_ft = 0;
		this.idx_tmp = 0;
		this.cap_vet = new CList[3];
		this.cup_vet = new CList[3];
		for (int i = 0; i < 3; i++) {
			this.cap_vet[i] = new CList();
			this.cap_vet[i].init(order.gene_num);
			this.cup_vet[i] = new CList();
			this.cup_vet[i].init(order.gene_num);
		}
		a = new CList();
		b = new CList();
		ap = new CList();
		bp = new CList();
		a.init(order.gene_num);
		b.init(order.gene_num);
		ap.init(order.gene_num);
		bp.init(order.gene_num);
		this.gene_order_to_graph(order);
		this.get_bounds();
	}

	public void init(Graph g) {
		// TODO Auto-generated constructor stub
		this.size = g.size;
		this.v_num = g.v_num;
		this.gene_num = g.gene_num;
		this.node_num = g.node_num;
		this.edge_num = g.edge_num;
		this.chr_num = g.chr_num;
		this.path00 = new int[3];
		this.path11 = new int[3];
		this.cycle00 = new int[3];
		this.cycle_ambi = new int[3];
		this.path_final = new int[3];
		for (int i = 0; i < 3; i++) {
			this.path00[i] = g.path00[i];
			this.path11[i] = g.path11[i];
			this.cycle00[i] = g.cycle00[i];
			this.cycle_ambi[i] = g.cycle_ambi[i];
			this.path_final[i] = g.path_final[i];
		}
		adj_mat = new int[v_num][3];
		check = new boolean[v_num];
		footprint = new int[gene_num * 4];
		major_tmp = new int[v_num];
		for (int i = 0; i < v_num; i++) {
			check[i] = g.check[i];
			footprint[i] = g.footprint[i];
			major_tmp[i] = g.major_tmp[i];
			for (int j = 0; j < 3; j++)
				adj_mat[i][j] = g.adj_mat[i][j];
		}

		pa_count = new int[v_num][3];
		vet_type = new int[v_num][3];
		vet_panum = new int[v_num][3];
		spec_path = new int[v_num][3];
		pa_path = new int[v_num][3];
		vet_rank = new int[v_num][3];
		vet_cynum = new int[v_num][3];

		this.cynum = new int[v_num][3];
		this.cyrank = new int[v_num][3];

		cycle = new int[3];
		cycle_number = g.cycle_number;
		c = new int[3];
		lower_bound = g.upper_bound;
		upper_bound = g.lower_bound;

		this.unused = new boolean[v_num];
		this.idx_ft = g.idx_ft;
		this.idx_tmp = g.idx_tmp;
		this.cap_vet = new CList[3];
		this.cup_vet = new CList[3];
		for (int i = 0; i < 3; i++) {
			this.cap_vet[i] = new CList();
			this.cap_vet[i].init(g.gene_num);
			this.cap_vet[i].init(g.cap_vet[i]);
			this.cup_vet[i] = new CList();
			this.cup_vet[i].init(g.gene_num);
			this.cup_vet[i].init(g.cup_vet[i]);
		}
		a = new CList();
		b = new CList();
		ap = new CList();
		bp = new CList();
		a.init(g.gene_num);
		b.init(g.gene_num);
		ap.init(g.gene_num);
		bp.init(g.gene_num);
	}

	@Override
	public void init(AdjNode n1, AdjNode n2, AdjNode n3) {
		// TODO Auto-generated constructor stub
		this.size = n1.adj.num_gene + n1.adj.num_chr;
		this.v_num = n1.adj.num_gene * 2;
		this.gene_num = n1.adj.num_gene;
		this.node_num = n1.adj.num_gene * 2;
		this.edge_num = n1.adj.num_gene;
		this.chr_num = n1.adj.num_chr;
		this.path00 = new int[3];
		this.path11 = new int[3];
		this.cycle00 = new int[3];
		this.cycle_ambi = new int[3];
		this.path_final = new int[3];
		adj_mat = new int[v_num][3];
		check = new boolean[v_num];
		footprint = new int[gene_num * 4];
		major_tmp = new int[v_num];
		for (int i = 0; i < v_num; i++) {
			check[i] = true;
			adj_mat[i][0] = n1.adj.reg_adj[i];
			adj_mat[i][1] = n2.adj.reg_adj[i];
			adj_mat[i][2] = n3.adj.reg_adj[i];
		}

		pa_count = new int[v_num][3];
		vet_type = new int[v_num][3];
		vet_panum = new int[v_num][3];
		spec_path = new int[v_num][3];
		pa_path = new int[v_num][3];
		vet_rank = new int[v_num][3];
		vet_cynum = new int[v_num][3];

		this.cynum = new int[v_num][3];
		this.cyrank = new int[v_num][3];

		cycle = new int[3];
		cycle_number = 0;
		c = new int[3];
		lower_bound = 0;
		upper_bound = 0;

		this.unused = new boolean[v_num];
		this.idx_ft = 0;
		this.idx_tmp = 0;
		this.cap_vet = new CList[3];
		this.cup_vet = new CList[3];
		// initialize cap nodes
		this.cap_vet[0] = new CList();
		this.cap_vet[0].init(gene_num);
		this.cap_vet[0].init(n1.adj.cap_adj);
		this.cup_vet[0] = new CList();
		this.cup_vet[0].init(gene_num);
		this.cap_vet[1] = new CList();
		this.cap_vet[1].init(gene_num);
		this.cap_vet[1].init(n2.adj.cap_adj);
		this.cup_vet[1] = new CList();
		this.cup_vet[1].init(gene_num);
		this.cap_vet[2] = new CList();
		this.cap_vet[2].init(gene_num);
		this.cap_vet[2].init(n3.adj.cap_adj);
		this.cup_vet[2] = new CList();
		this.cup_vet[2].init(gene_num);

		a = new CList();
		b = new CList();
		ap = new CList();
		bp = new CList();
		a.init(gene_num);
		b.init(gene_num);
		ap.init(gene_num);
		bp.init(gene_num);

		// initialize the adj
		ft_before_idx = 0;
		map = new int[gene_num * 2];
		map_o = new int[gene_num * 2];
		for (int i = 0; i < gene_num * 2; i++) {
			map[i]=i;
			map_o[i]=i;
		}
		ft_before_rename = new int[gene_num * 2];
		median_adj = new Adjacency(gene_num, 1, "");
		for (int i = 0; i < gene_num * 2; i++) {
			median_adj.reg_adj[i] = Constant.NULL;
		}
		this.get_bounds();
	}

	@Override
	public void gene_order_to_graph(GeneOrder o) {
		// TODO Auto-generated method stub
		int color, gene;
		int i;
		int order, order_forward, order_backward;
		int forward_pos, backward_pos;
		int forward_pt, backward_pt;

		for (color = 0; color < 3; color++) {
			for (gene = 1; gene < (gene_num + o.chr_num[color]); gene++) {
				// (forward_pt)<-(forward_pos) (backward_pos)->(backward_pt)
				order = o.gene_order[gene][color];
				if ((Math.abs(order)) != Const.CAP) {
					// forward
					forward_pos = order > 0 ? (Math.abs(order) - 1) * 2 : (Math
							.abs(order) - 1) * 2 + 1;
					order_forward = o.gene_order[gene - 1][color];
					if (Math.abs(order_forward) != Const.CAP) {
						forward_pt = order_forward > 0 ? (Math
								.abs(order_forward) - 1) * 2 + 1 : (Math
								.abs(order_forward) - 1) * 2;
						adj_mat[forward_pos][color] = forward_pt;
					} else {
						adj_mat[forward_pos][color] = Const.CAP;
						cap_vet[color].add(forward_pos);
					}

					// backward
					backward_pos = order > 0 ? (Math.abs(order) - 1) * 2 + 1
							: (Math.abs(order) - 1) * 2;
					order_backward = o.gene_order[gene + 1][color];
					if (Math.abs(order_backward) != Const.CAP) {
						backward_pt = order_backward > 0 ? (Math
								.abs(order_backward) - 1) * 2 : (Math
								.abs(order_backward) - 1) * 2 + 1;
						adj_mat[backward_pos][color] = backward_pt;
					} else {
						adj_mat[backward_pos][color] = Const.CAP;
						cap_vet[color].add(backward_pos);
					}

				} else if (o.gene_order[gene - 1][color] == Const.CAP
						&& gene == 1)
					path00[color]++;
				else if (o.gene_order[gene + 1][color] == Const.CAP)
					path00[color]++;
			}
		}

		cycle_number = 0;
		size = gene_num + o.chr_num[0];
		for (i = 0; i < node_num; i++)
			check[i] = true;
		// print_adj_mat();
	}

	@Override
	public void shrink(int[] major, int start, int end) {
		// TODO Auto-generated method stub
		int m = (end - start) / 2;
		int cap_edges = 0;
		for (int i = start; i < end; i += 2) {
			int l = major[i]; // the left vertex
			int r = major[i + 1]; // the right vertex

			if (l != Const.CAP && r != Const.CAP) { // not connected to CAP
				check[l] = false;
				check[r] = false;
				for (int color = 0; color < 3; color++) {
					if (adj_mat[l][color] == r) {
						cycle_number++;
					} else {
						int lc = adj_mat[l][color];
						int rc = adj_mat[r][color];
						if (lc != Const.CAP && lc != Const.CUP
								&& rc != Const.CAP && rc != Const.CUP) {
							adj_mat[lc][color] = rc;
							adj_mat[rc][color] = lc;
						} else if (lc == Const.CAP && rc != Const.CAP
								&& rc != Const.CUP) {
							adj_mat[rc][color] = Const.CAP;
							cap_vet[color].remove(l);
							cap_vet[color].add(rc);
						} else if (lc == Const.CUP && rc != Const.CAP
								&& rc != Const.CUP) {
							adj_mat[rc][color] = Const.CUP;
							cup_vet[color].remove(l);
							cup_vet[color].add(rc);
							if (rc >= Const.CUP)
								System.out.printf("sh1");
						} else if (lc != Const.CAP && lc != Const.CUP
								&& rc == Const.CAP) {
							adj_mat[lc][color] = Const.CAP;
							cap_vet[color].remove(r);
							cap_vet[color].add(lc);
						} else if (lc != Const.CAP && lc != Const.CUP
								&& rc == Const.CUP) {
							adj_mat[lc][color] = Const.CUP;
							cup_vet[color].remove(r);
							cup_vet[color].add(lc);
						} else if (lc == Const.CAP && rc == Const.CAP) {
							cap_vet[color].remove(l);
							cap_vet[color].remove(r);
							if (path11[color] > 0) {
								cycle_number++;
								cycle_ambi[color]++;
								path11[color]--;
							} else
								path00[color]++;
						} else if (lc == Const.CUP && rc == Const.CUP) {
							cup_vet[color].remove(l);
							cup_vet[color].remove(r);
							if (path00[color] > 0) {
								cycle_number++;
								cycle_ambi[color]++;
								path00[color]--;
							} else
								path11[color]++;
						} else if (lc == Const.CUP && rc == Const.CAP) {
							cup_vet[color].remove(l);
							cap_vet[color].remove(r);
							cycle_number++;
						} else if (lc == Const.CAP && rc == Const.CUP) {
							cap_vet[color].remove(l);
							cup_vet[color].remove(r);
							cycle_number++;
						}
					}
				}

			} else if (l != Const.CAP && r == Const.CAP) {
				check[l] = false;
				for (int color = 0; color < 3; color++) {
					int lc = adj_mat[l][color];
					if (lc == Const.CAP) {
						cycle_number++;
						cap_vet[color].remove(l);
					} else if (lc == Const.CUP) {
						cup_vet[color].remove(l);
						if (path00[color] > 0) {
							path00[color]--;
							cycle_ambi[color]++;
							cycle_number++;
						} else
							path11[color]++;
					} else {
						adj_mat[lc][color] = Const.CUP;
						cup_vet[color].add(lc);
					}
				}
				cap_edges++;

			} else if (l == Const.CAP && r == Const.CAP) {
				for (int color = 0; color < 3; color++) {
					if (path00[color] > 0) {
						path00[color]--;
						cycle_ambi[color]++;
						cycle_number++;
					} else
						path11[color]++;
				}
				cap_edges += 2;
			} else {
				System.out.printf("error left=Const.CAP and right!=Const.CAP");
				System.exit(1);
			}
		}

		size -= m;
		node_num = node_num - (end - start) + cap_edges;
		chr_num -= cap_edges;
		// if (chr_num < 0) {
		// System.out.printf("chr_num is smaller than 0");
		// exit(1);
		// }

		if (chr_num > 0 && node_num == 0) {
			chr_num_final = chr_num;
			for (int color = 0; color < 3; color++) {
				if (chr_num / 2 != size) {
					System.out.printf("size is note equal to chr_num/2");
					//System.exit(1);
				}
				path11[color] += chr_num / 2;
				if (path00[color] > 0) {
					int min = path00[color] > path11[color] ? path11[color]
							: path00[color];
					if (path00[color] != path11[color]) {
						System.out.printf("different number of paths");
						//System.exit(1);
					}
					cycle_number += min;
					cycle_ambi[color] += min;
					path00[color] -= min;
					path11[color] -= min;
					path_final[color] = min;
				}
			}
			size = chr_num = node_num = 0;
		}
		// updates footprints
		for (int j = start; j < end; j++)
			this.add_ft(major[j]);
	}

	public void expand(int[] major, int start, int end) {
		// TODO Auto-generated method stub
		int m = (end - start) / 2; // the number of black edges, i.e. the size
		// of the subgraph

		// resume the previous step;
		if (chr_num_final > 0 && node_num == 0) {
			for (int color = 0; color < 3; color++) {
				path11[color] -= chr_num_final / 2;
				if (path_final[color] > 0) {
					cycle_number -= path_final[color];
					cycle_ambi[color] -= path_final[color];
					path00[color] += path_final[color];
					path11[color] += path_final[color];
				}
			}
			node_num = 0;
			size = chr_num / 2;
			chr_num = chr_num_final;
		}

		// refactor the adj_mat array and the cap array
		int cap_edges = 0;
		for (int i = end - 1; i >= start; i -= 2) {
			int l = major[i - 1]; // the left vertex
			int r = major[i]; // the right vertex

			if (l != Const.CAP && r != Const.CAP) { // not connected to
													// Const.CAP
				check[l] = true;
				check[r] = true;
				for (int color = 0; color < 3; color++) {
					if (adj_mat[l][color] == r) {
						cycle_number--;
					} else {
						int lc = adj_mat[l][color];
						int rc = adj_mat[r][color];
						if (lc != Const.CAP && lc != Const.CUP
								&& rc != Const.CAP && rc != Const.CUP) {
							adj_mat[lc][color] = l;
							adj_mat[rc][color] = r;
						} else if (lc == Const.CAP && rc != Const.CAP
								&& rc != Const.CUP) {
							adj_mat[rc][color] = r;
							cap_vet[color].remove(rc);
							cap_vet[color].add(l);
						} else if (lc == Const.CUP && rc != Const.CAP
								&& rc != Const.CUP) {
							adj_mat[rc][color] = r;
							cup_vet[color].remove(rc);
							cup_vet[color].add(l);
						} else if (lc != Const.CAP && lc != Const.CUP
								&& rc == Const.CAP) {
							adj_mat[lc][color] = l;
							cap_vet[color].remove(lc);
							cap_vet[color].add(r);
						} else if (lc != Const.CAP && lc != Const.CUP
								&& rc == Const.CUP) {
							adj_mat[lc][color] = l;
							cup_vet[color].remove(lc);
							cup_vet[color].add(r);
						} else if (lc == Const.CAP && rc == Const.CAP) {
							cap_vet[color].add(l);
							cap_vet[color].add(r);
							if (cycle_ambi[color] > 0) {
								cycle_number--;
								path11[color]++;
								cycle_ambi[color]--;
							} else
								path00[color]--;
						} else if (lc == Const.CUP && rc == Const.CUP) {
							cup_vet[color].add(l);
							cup_vet[color].add(r);
							if (cycle_ambi[color] > 0) {
								cycle_number--;
								cycle_ambi[color]--;
								path00[color]++;
							} else
								path11[color]--;
						} else if (lc == Const.CUP && rc == Const.CAP) {
							cup_vet[color].add(l);
							cap_vet[color].add(r);
							cycle_number--;
						} else if (lc == Const.CAP && rc == Const.CUP) {
							cap_vet[color].add(l);
							cup_vet[color].add(r);
							cycle_number--;
						}
					}
				}

			} else if (l != Const.CAP && r == Const.CAP) {
				check[l] = true;
				for (int color = 0; color < 3; color++) {
					int lc = adj_mat[l][color];
					if (lc == Const.CAP) {
						cycle_number--;
						cap_vet[color].add(l);
					} else if (lc == Const.CUP) {
						cup_vet[color].add(l);
						if (cycle_ambi[color] > 0) {
							path00[color]++;
							cycle_number--;
							cycle_ambi[color]--;
						} else
							path11[color]--;
					} else {
						adj_mat[lc][color] = l;
						cup_vet[color].remove(lc);
					}
				}
				cap_edges--;

			} else if (l == Const.CAP && r == Const.CAP) {
				for (int color = 0; color < 3; color++) {
					if (cycle00[color] > 0) {
						path00[color]++;
						cycle_number--;
						cycle00[color]--;
					} else
						path11[color]--;
				}
				cap_edges -= 2;
			} else {
				System.exit(1);
			}
		}
		size += m;
		node_num = node_num + (end - start) + cap_edges;
		chr_num -= cap_edges;
		// updates footprints
		for (int j = start; j < end; j++)
			this.mv_back_ft();
	}

	@Override
	public void rename() {
		// TODO Auto-generated method stub
		cap_vet[0].clear();
		cap_vet[1].clear();
		cap_vet[2].clear();
		cup_vet[0].clear();
		cup_vet[1].clear();
		cup_vet[2].clear();

		int new_num_vet = 0;
		int idx = 0;
		map = new int[Const.MAX_GENE_NUM * 2];
		map_o = new int[Const.MAX_GENE_NUM * 2];
		for (int i = 0; i < v_num; i++) {
			if (check[i] == true) {
				map[i] = idx;
				map_o[idx] = i;
				idx++;
				new_num_vet++;
			} else
				map[i] = -1;
		}
		idx = 0;
		for (int i = 0; i < v_num; i++) {
			if (check[i] == true) {
				int c1 = adj_mat[i][0];
				int c2 = adj_mat[i][1];
				int c3 = adj_mat[i][2];
				if (c1 < Const.CUP)
					adj_mat[idx][0] = map[c1];
				else if (c1 == Const.CAP) {
					adj_mat[idx][0] = c1;
					cap_vet[0].add(idx);
				} else if (c1 == Const.CUP) {
					adj_mat[idx][0] = c1;
					cup_vet[0].add(idx);
				}

				if (c2 < Const.CUP)
					adj_mat[idx][1] = map[c2];
				else if (c2 == Const.CAP) {
					adj_mat[idx][1] = c2;
					cap_vet[1].add(idx);
				} else if (c2 == Const.CUP) {
					adj_mat[idx][1] = c2;
					cup_vet[1].add(idx);
				}

				if (c3 < Const.CUP)
					adj_mat[idx][2] = map[c3];
				else if (c3 == Const.CAP) {
					adj_mat[idx][2] = c3;
					cap_vet[2].add(idx);
				} else if (c3 == Const.CUP) {
					adj_mat[idx][2] = c3;
					cup_vet[2].add(idx);
				}
				idx++;
			}
		}

		for (int i = 0; i < new_num_vet; i++)
			check[i] = true;

		node_num = new_num_vet;
		v_num = new_num_vet;
		gene_num = new_num_vet / 2;
		edge_num = new_num_vet / 2;
	}

	@Override
	public Adjacency toMedianAdj() {
		// copy the renamed vertices
		for (int i = 0; i < ft_before_idx; i += 2) {
			if (ft_before_rename[i] < Const.CUP
					&& ft_before_rename[i + 1] < Const.CUP) {
				median_adj.reg_adj[ft_before_rename[i]] = ft_before_rename[i + 1];
				median_adj.reg_adj[ft_before_rename[i + 1]] = ft_before_rename[i];
			} else if (ft_before_rename[i] >= Const.CUP) {
				median_adj.reg_adj[ft_before_rename[i + 1]] = ft_before_rename[i];
				median_adj.cap_adj.add(ft_before_rename[i + 1]);
			} else if (ft_before_rename[i + 1] >= Const.CUP) {
				median_adj.reg_adj[ft_before_rename[i]] = ft_before_rename[i + 1];
				median_adj.cap_adj.add(ft_before_rename[i]);
			}
		}
		// copy the current ft vertices
		for (int i = 0; i < this.idx_ft; i += 2) {
			if(footprint[i]>=Const.CUP || footprint[i+1]>=Const.CUP)
				continue;
			if (map_o[footprint[i]] < Const.CUP
					&& map_o[footprint[i + 1]] < Const.CUP) {
				median_adj.reg_adj[map_o[footprint[i]]] = map_o[footprint[i + 1]];
				median_adj.reg_adj[map_o[footprint[i + 1]]] = map_o[footprint[i]];
			} else if (map_o[footprint[i]] >= Const.CUP) {
				median_adj.reg_adj[map_o[footprint[i + 1]]] = map_o[footprint[i]];
				median_adj.cap_adj.add(map_o[footprint[i + 1]]);
			} else if (map_o[footprint[i + 1]] >= Const.CUP) {
				median_adj.reg_adj[map_o[footprint[i]]] = map_o[footprint[i + 1]];
				median_adj.cap_adj.add(map_o[footprint[i]]);
			}
		}

		int from_adj = -1;
		int to_adj = -1;
		// greedy add cap
		if (median_adj.cap_adj.size() == 0) {
			from_adj = this.next_median_adj();
			to_adj = this.next_median_adj();
			System.out.println(from_adj + " " + to_adj);
			median_adj.reg_adj[from_adj] = Const.CAP;
			median_adj.cap_adj.add(from_adj);
			median_adj.reg_adj[to_adj] = Const.CAP;
			median_adj.cap_adj.add(to_adj);
		} else if (median_adj.cap_adj.size() == 1) {
			from_adj = this.next_median_adj();
			median_adj.reg_adj[from_adj] = Const.CAP;
			median_adj.cap_adj.add(from_adj);
		}
		// greedy add nodes
		while ((from_adj = this.next_median_adj()) != -1) {
			to_adj = this.next_median_adj();
			median_adj.reg_adj[from_adj] = to_adj;
			median_adj.reg_adj[to_adj] = from_adj;
		}
		return this.median_adj;
	}

	@Override
	public void get_bounds() {
		// TODO Auto-generated method stub
		count_cycle(0, 1);
		count_cycle(0, 2);
		count_cycle(1, 2);

		upper_bound = cycle_number
				+ (int) Math.floor((3 * size + c[2] + c[1] + c[0]) / 2.0);
		int tc[] = { c[2] + c[1], c[2] + c[0], c[1] + c[0] };
		int tcm = tc[0], index = 0;
		for (int i = 1; i < 3; i++) {
			if (tc[i] > tcm) {
				tcm = tc[i];
				index = i;
			}
		}
		lower_bound = cycle_number + tcm;
	}

	@Override
	public void count_cycle(int c1, int c2) {
		// TODO Auto-generated method stub
		int i, cycles;
		int start, left, right;
		int co;

		int pathaa = path00[c1]; //
		int pathbb = path00[c2]; //
		int pathab = 0; //
		int pathapap = path11[c1];
		int pathbpbp = path11[c2];
		int pathapbp = 0;
		int pathaap = 0; //
		int pathbbp = 0; //
		int pathabp = 0; //
		int pathbap = 0; //

		a.init(cap_vet[c1]);
		b.init(cap_vet[c2]);
		ap.init(cup_vet[c1]);
		bp.init(cup_vet[c2]);

		if (c1 == 0 && c2 == 1)
			co = 2;
		else if (c1 == 0 && c2 == 2)
			co = 1;
		else
			co = 0;

		cycles = 0;
		for (i = 0; i < v_num; i++) {
			if (check[i] == true)
				unused[i] = true;
			else
				unused[i] = false;
		}

		// this is to calculate path number
		while (a.size() != 0) {
			left = a.pop_back();
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c2];
				if (right == Const.CAP) {
					pathab++;
					b.remove(left);
					break;
				} else if (right == Const.CUP) {
					pathabp++;
					bp.remove(left);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c1];
				if (left == Const.CAP) {
					pathaa++;
					a.remove(right);
					break;
				} else if (left == Const.CUP) {
					pathaap++;
					ap.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		// this is also to calculate path number
		while (b.size() != 0) {
			left = b.pop_back();
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c1];
				if (right == Const.CUP) {
					pathbap++;
					ap.remove(left);
					break;
				}
				if (right == Const.CAP) {
					pathaa++;
					a.remove(left);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c2];
				if (left == Const.CAP) {
					pathbb++;
					b.remove(right);
					break;
				} else if (left == Const.CUP) {
					pathbbp++;
					bp.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		// this is the third path calculator
		while (ap.size() != 0) {
			left = ap.pop_back();
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c2];
				if (right == Const.CUP) {
					pathapbp++;
					bp.remove(left);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c1];
				if (left == Const.CUP) {
					pathapap++;
					ap.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		// this is the fourth path calculator
		while (bp.size() != 0) {
			left = bp.pop_back();
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c1];
				unused[right] = false;

				left = adj_mat[right][c2];
				if (left == Const.CUP) {
					pathbpbp++;
					bp.remove(right);
					break;
				}
				if (left == Const.CAP) {
					pathbpbp++;
					a.remove(right);
					break;
				}
				unused[left] = false;
			}
		}

		// this is to calculate cycle numbers
		for (i = 0; i < v_num; i++) {
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
		pathapbp += 2 * chr_num;
		c[co] = cycles + pathaap + pathbbp
				+ (pathaa < pathapap ? pathaa : pathapap)
				+ (pathbb < pathbpbp ? pathbb : pathbpbp)
				+ (pathabp + pathbap + pathab + pathapbp) / 2;
	}

	@Override
	public void get_rank_cynum(int c1, int c2) {
		// TODO Auto-generated method stub
		int i, cycles;
		int start, left, right;
		int pa_num = 0;

		int co;
		int rank;

		cycles = 0;

		int pathaa = 0; //
		int pathbb = 0; //
		int pathab = 0; //
		int pathapap = 0;
		int pathbpbp = 0;
		int pathapbp = 0;
		int pathaap = 0; //
		int pathbbp = 0; //
		int pathabp = 0; //
		int pathbap = 0; //

		if (c1 == 0 && c2 == 1)
			co = 2;
		else if (c1 == 0 && c2 == 2)
			co = 1;
		else
			co = 0;

		for (i = 0; i < v_num; i++) {
			if (check[i] == true)
				unused[i] = true;
			else
				unused[i] = false;
		}

		a.init(cap_vet[c1]);
		b.init(cap_vet[c2]);
		ap.init(cup_vet[c1]);
		bp.init(cup_vet[c2]);

		// this is to calculate path number
		while (a.size() != 0) {
			pa_count[pa_num][co] = 0;
			left = a.pop_back();
			set_vet(left, co, Const.TYPE_PATH, pa_num);
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c2];
				if (right < Const.CUP)
					set_vet(right, co, Const.TYPE_PATH, pa_num);
				if (right == Const.CAP) {
					pathab++;
					b.remove(left);
					pa_num = set_path(pa_num, co, Const.AB);
					break;
				} else if (right == Const.CUP) {
					pathabp++;
					bp.remove(left);
					pa_num = set_path(pa_num, co, Const.ABP);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c1];
				if (left < Const.CUP)
					set_vet(left, co, Const.TYPE_PATH, pa_num);
				if (left == Const.CAP) {
					pathaa++;
					a.remove(right);
					pa_num = set_path(pa_num, co, Const.AA);
					break;
				} else if (left == Const.CUP) {
					pathaap++;
					ap.remove(right);
					pa_num = set_path(pa_num, co, Const.AAP);
					break;
				}
				unused[left] = false;
			}
		}
		// this is also to calculate path number
		while (b.size() != 0) {
			pa_count[pa_num][co] = 0;
			left = b.pop_back();
			set_vet(left, co, Const.TYPE_PATH, pa_num);
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c1];
				if (right < Const.CUP)
					set_vet(right, co, Const.TYPE_PATH, pa_num);
				if (right == Const.CUP) {
					pathbap++;
					ap.remove(left);
					pa_num = set_path(pa_num, co, Const.BAP);
					break;
				}

				unused[right] = false;

				left = adj_mat[right][c2];
				if (left < Const.CUP)
					set_vet(left, co, Const.TYPE_PATH, pa_num);

				if (left == Const.CAP) {
					pathbb++;
					b.remove(right);
					pa_num = set_path(pa_num, co, Const.BB);
					break;
				} else if (left == Const.CUP) {
					pathbbp++;
					bp.remove(right);
					pa_num = set_path(pa_num, co, Const.BBP);
					break;
				}
				unused[left] = false;
			}
		}
		// this is the third path calculator
		while (ap.size() != 0) {
			pa_count[pa_num][co] = 0;
			left = ap.pop_back();
			set_vet(left, co, Const.TYPE_PATH, pa_num);
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c2];
				if (right < Const.CUP)
					set_vet(right, co, Const.TYPE_PATH, pa_num);
				if (right == Const.CUP) {
					pathapbp++;
					bp.remove(left);
					pa_num = set_path(pa_num, co, Const.APBP);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c1];
				if (left < Const.CUP)
					set_vet(left, co, Const.TYPE_PATH, pa_num);

				if (left == Const.CUP) {
					pathapap++;
					ap.remove(right);
					pa_num = set_path(pa_num, co, Const.APAP);
					break;
				}
				unused[left] = false;
			}
		}
		// this is the fourth path calculator
		while (bp.size() != 0) {
			pa_count[pa_num][co] = 0;
			left = bp.pop_back();
			set_vet(left, co, Const.TYPE_PATH, pa_num);
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c1];
				if (right < Const.CUP)
					set_vet(right, co, Const.TYPE_PATH, pa_num);
				unused[right] = false;

				left = adj_mat[right][c2];
				if (left < Const.CUP)
					set_vet(left, co, Const.TYPE_PATH, pa_num);
				if (left == Const.CUP) {
					pathbpbp++;
					bp.remove(right);
					pa_num = set_path(pa_num, co, Const.BPBP);
					break;
				}
				unused[left] = false;
			}
		}
		// this is to calculate cycles
		for (i = 0; i < v_num; i++) {
			if (unused[i] != true)
				continue;
			start = left = i;
			rank = 0;
			do {
				right = adj_mat[left][c1];
				vet_type[right][co] = 1;
				unused[left] = unused[right] = false;
				vet_cynum[left][co] = cycles;
				vet_panum[left][co] = pa_num;
				pa_count[pa_num][co] += 1;
				vet_rank[left][co] = rank;
				vet_cynum[right][co] = cycles;
				vet_panum[right][co] = pa_num;
				pa_count[pa_num][co] += 1;
				vet_rank[right][co] = rank + 1;
				rank += 2;
				left = adj_mat[right][c2];
				vet_type[left][co] = 1;
			} while (left != start);
			cycles++;
			pa_num++;
		}

		// pathapbp+=chr_num;

		spec_path[Const.CYCLE][co] = cycles;
		spec_path[Const.AAP][co] = pathaap;
		spec_path[Const.BBP][co] = pathbbp;
		spec_path[Const.ABP][co] = pathabp;
		spec_path[Const.BAP][co] = pathbap;
		spec_path[Const.AB][co] = pathab;
		spec_path[Const.APBP][co] = pathapbp;
		spec_path[Const.AA][co] = pathaa;
		spec_path[Const.BB][co] = pathbb;
		spec_path[Const.APAP][co] = pathapap;
		spec_path[Const.BPBP][co] = pathbpbp;
	}

	void set_vet(int vet_id, int co, int type, int panum) {
		vet_type[vet_id][co] = type;
		vet_panum[vet_id][co] = panum;
		pa_count[panum][co] += 1;
	}

	int set_path(int pa_num, int co, int type) {
		this.pa_path[pa_num][co] = type;
		return ++pa_num;
	}

	@Override
	public void get_bounds_linear(int v1, int v2) {
		// TODO Auto-generated method stub
		count_cycle_linear(v1, v2, 0, 1);
		count_cycle_linear(v1, v2, 0, 2);
		count_cycle_linear(v1, v2, 1, 2);

		upper_bound = cycle_number
				+ (int) Math.floor((3 * size + c[2] + c[1] + c[0]) / 2.0);
		int tc[] = { c[2] + c[1], c[2] + c[0], c[1] + c[0] };
		int tcm = tc[0], index = 0;
		for (int i = 1; i < 3; i++) {
			if (tc[i] > tcm) {
				tcm = tc[i];
				index = i;
			}
		}
		lower_bound = cycle_number + tcm;
	}

	@Override
	public void count_cycle_linear(int v1, int v2, int c1, int c2) {
		// TODO Auto-generated method stub
		int color;
		int single_path[][] = { { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 },
				{ 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 }, { 0, 0 },
				{ 0, 0 } };
		int left, right;
		int type1, type2;
		int comp1, comp2;

		if (c1 == 0 && c2 == 1)
			color = 2;
		else if (c1 == 0 && c2 == 2)
			color = 1;
		else
			color = 0;

		int cycle = spec_path[Const.CYCLE][color];
		int sum = 0;

		if (v2 >= Const.CUP) {
			type1 = vet_type[v1][color];
			comp1 = vet_panum[v1][color];
			type2 = -1;
			comp2 = -1;
			if (type1 == Const.TYPE_CYCLE) {
				cycle--;
			} else {
				single_path[pa_path[comp1][color]][0] = 1;
			}
		} else {
			type1 = vet_type[v1][color];
			type2 = vet_type[v2][color];
			comp1 = vet_panum[v1][color];
			comp2 = vet_panum[v2][color];

			if (type1 == Const.TYPE_CYCLE && type2 == Const.TYPE_CYCLE) {
				update_on_cycles(v1, v2, color, c1, c2);
				return;
			} else if (type1 == Const.TYPE_CYCLE && type2 != Const.TYPE_CYCLE) {
				single_path[pa_path[comp2][color]][0]++;
				cycle--;
			} else if (type1 != Const.TYPE_CYCLE && type2 == Const.TYPE_CYCLE) {
				single_path[pa_path[comp1][color]][0]++;
				cycle--;
			} else if (type1 != Const.TYPE_CYCLE && type2 != Const.TYPE_CYCLE
					&& comp1 == comp2) {
				single_path[pa_path[comp1][color]][0]++;
				sum = pa_count[comp1][color] - 2;
			} else if (type1 != Const.TYPE_CYCLE && type2 != Const.TYPE_CYCLE
					&& comp1 != comp2) {
				single_path[pa_path[comp1][color]][0]++;
				single_path[pa_path[comp2][color]][0]++;
				sum = pa_count[comp2][color] + pa_count[comp1][color] - 2;
			}
		}

		for (int i = 0; i < cap_vet[c1].size(); i++) {
			if (vet_panum[cap_vet[c1].c[i]][color] == comp1
					|| vet_panum[cap_vet[c1].c[i]][color] == comp2)
				a.add(cap_vet[c1].c[i]);
		}
		for (int i = 0; i < cap_vet[c2].size(); i++) {
			if (vet_panum[cap_vet[c2].c[i]][color] == comp1
					|| vet_panum[cap_vet[c2].c[i]][color] == comp2)
				b.add(cap_vet[c2].c[i]);
		}
		for (int i = 0; i < cup_vet[c1].size(); i++) {
			if (vet_panum[cup_vet[c1].c[i]][color] == comp1
					|| vet_panum[cup_vet[c1].c[i]][color] == comp2)
				ap.add(cup_vet[c1].c[i]);
		}
		for (int i = 0; i < cup_vet[c2].size(); i++) {
			if (vet_panum[cup_vet[c2].c[i]][color] == comp1
					|| vet_panum[cup_vet[c2].c[i]][color] == comp2)
				bp.add(cup_vet[c2].c[i]);
		}

		while (a.size() != 0) {
			left = a.pop_back();
			sum--;
			while (true) {
				right = adj_mat[left][c2];
				if (right < Const.CUP)
					sum--;
				if (right == Const.CAP) {
					single_path[Const.AB][1]++;
					b.remove(left);
					break;
				} else if (right == Const.CUP) {
					single_path[Const.ABP][1]++;
					bp.remove(left);
					break;
				}

				left = adj_mat[right][c1];
				if (left < Const.CUP)
					sum--;
				if (left == Const.CAP) {
					single_path[Const.AA][1]++;
					a.remove(right);
					break;
				} else if (left == Const.CUP) {
					single_path[Const.AAP][1]++;
					ap.remove(right);
					break;
				}
			}
		}
		// this is also to calculate path number
		while (b.size() != 0) {
			left = b.pop_back();
			sum--;
			while (true) {
				right = adj_mat[left][c1];
				if (right < Const.CUP)
					sum--;
				if (right == Const.CUP) {
					single_path[Const.BAP][1]++;
					ap.remove(left);
					break;
				}

				left = adj_mat[right][c2];
				if (left < Const.CUP)
					sum--;
				if (left == Const.CAP) {
					single_path[Const.BB][1]++;
					b.remove(right);
					break;
				} else if (left == Const.CUP) {
					single_path[Const.BBP][1]++;
					bp.remove(right);
					break;
				}
			}
		}
		// this is the third path calculator
		while (ap.size() != 0) {
			left = ap.pop_back();
			sum--;
			while (true) {
				right = adj_mat[left][c2];
				if (right < Const.CUP)
					sum--;
				if (right == Const.CUP) {
					single_path[Const.APBP][1]++;
					bp.remove(left);
					break;
				}

				left = adj_mat[right][c1];
				if (left < Const.CUP)
					sum--;

				if (left == Const.CUP) {
					single_path[Const.APAP][1]++;
					ap.remove(right);
					break;
				}
			}
		}
		// this is the fourth path calculator
		while (bp.size() != 0) {
			left = bp.pop_back();
			sum--;
			while (true) {
				right = adj_mat[left][c1];
				if (right < Const.CUP)
					sum--;

				left = adj_mat[right][c2];
				if (left < Const.CUP)
					sum--;
				if (left == Const.CUP) {
					single_path[Const.BPBP][1]++;
					bp.remove(right);
					break;
				}
			}
		}

		if (sum > 0)
			cycle++;

		int aap = spec_path[Const.AAP][color] - single_path[Const.AAP][0]
				+ single_path[Const.AAP][1];
		int bbp = spec_path[Const.BBP][color] - single_path[Const.BBP][0]
				+ single_path[Const.BBP][1];
		int abp = spec_path[Const.ABP][color] - single_path[Const.ABP][0]
				+ single_path[Const.ABP][1];
		int bap = spec_path[Const.BAP][color] - single_path[Const.BAP][0]
				+ single_path[Const.BAP][1];
		int ab = spec_path[Const.AB][color] - single_path[Const.AB][0]
				+ single_path[Const.AB][1];
		int apbp = 2 * chr_num + spec_path[Const.APBP][color]
				- single_path[Const.APBP][0] + single_path[Const.APBP][1];
		int aa = spec_path[Const.AA][color] - single_path[Const.AA][0]
				+ single_path[Const.AA][1] + path00[c1];
		int bb = spec_path[Const.BB][color] - single_path[Const.BB][0]
				+ single_path[Const.BB][1] + path00[c2];
		int apap = spec_path[Const.APAP][color] - single_path[Const.APAP][0]
				+ single_path[Const.APAP][1] + path11[c1];
		int bpbp = spec_path[Const.BPBP][color] - single_path[Const.BPBP][0]
				+ single_path[Const.BPBP][1] + path11[c2];
		c[color] = cycle + aap + bbp + (abp + bap + apbp + ab) / 2
				+ (aa < apap ? aa : apap) + (bb < bpbp ? bb : bpbp);
	}

	void update_on_cycles(int v1, int v2, int color, int c1, int c2) {
		int x, y, w, z;
		int r_x, r_y, r_w, r_z;

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

	@Override
	public void get_start(int c1, int c2) {

		int i, cycles;
		int start, left, right;

		int rank;
		int color;

		if (c1 == 0 && c2 == 1)
			color = 2;
		else if (c1 == 0 && c2 == 2)
			color = 1;
		else
			color = 0;

		cycles = 0;

		for (i = 0; i < v_num; i++) {

			unused[i] = check[i];
		}

		a.init(cap_vet[c1]);
		b.init(cap_vet[c2]);
		ap.init(cup_vet[c1]);
		bp.init(cup_vet[c2]);

		// this is to calculate path number
		rank = 0;
		while (a.size() != 0) {
			left = a.pop_back();
			cynum[left][color] = cycles;
			rank++;
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c2];
				if (right < Const.CUP) {
					cynum[right][color] = cycles;
					rank++;
				}
				if (right == Const.CAP) {
					b.remove(left);
					break;
				} else if (right == Const.CUP) {
					bp.remove(left);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c1];
				if (left < Const.CUP) {
					cynum[left][color] = cycles;
					rank++;
				}
				if (left == Const.CAP) {
					a.remove(right);
					break;
				} else if (left == Const.CUP) {
					ap.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		if (rank > 0) {
			cyrank[cycles][color] = rank;
			cycles++;
		}

		// this is also to calculate path number
		rank = 0;
		while (b.size() != 0) {
			left = b.pop_back();
			cynum[left][color] = cycles;
			rank++;
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c1];
				if (right < Const.CUP) {
					cynum[right][color] = cycles;
					rank++;
				}
				if (right == Const.CUP) {
					ap.remove(left);
					break;
				}
				if (right == Const.CAP) {
					a.remove(left);
					b.remove(left);
					break;
				}

				unused[right] = false;

				left = adj_mat[right][c2];
				if (left < Const.CUP) {
					cynum[left][color] = cycles;
					rank++;
				}

				if (left == Const.CAP) {
					b.remove(right);
					break;
				} else if (left == Const.CUP) {
					bp.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		if (rank > 0) {
			cyrank[cycles][color] = rank;
			cycles++;
		}
		// this is the third path calculator
		rank = 0;
		while (ap.size() != 0) {
			left = ap.pop_back();
			cynum[left][color] = cycles;
			rank++;
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c2];
				if (right < Const.CUP) {
					cynum[right][color] = cycles;
					rank++;
				}
				if (right == Const.CUP) {
					bp.remove(left);
					break;
				}
				unused[right] = false;

				left = adj_mat[right][c1];
				if (left < Const.CUP) {
					cynum[left][color] = cycles;
					rank++;
				}

				if (left == Const.CUP) {
					ap.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		if (rank > 0) {
			cyrank[cycles][color] = rank;
			cycles++;
		}
		// this is the fourth path calculator
		rank = 0;
		while (bp.size() != 0) {
			left = bp.pop_back();
			cynum[left][color] = cycles;
			rank++;
			unused[left] = false;
			while (true) {
				right = adj_mat[left][c1];
				if (right < Const.CUP) {
					cynum[right][color] = cycles;
					rank++;
				}
				unused[right] = false;

				left = adj_mat[right][c2];
				if (left < Const.CUP) {
					cynum[left][color] = cycles;
					rank++;
				}
				if (left == Const.CUP) {
					bp.remove(right);
					break;
				}
				unused[left] = false;
			}
		}
		if (rank > 0) {
			cyrank[cycles][color] = rank;
			cycles++;
		}
		// cycle
		for (i = 0; i < v_num; i++) {
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

	public String toHash() {
		StringBuffer buffer = new StringBuffer();
		for (int i = 0; i < this.idx_ft; i++)
			buffer.append(this.footprint[i]);
		return buffer.toString();
	}

	public String toString() {
		StringBuffer buffer = new StringBuffer();
		int new_num_vet = 0;
		int idx = 0;
		int map[] = new int[v_num];
		for (int i = 0; i < v_num; i++) {
			if (check[i] == true) {
				map[i] = idx;
				idx++;
				new_num_vet++;
			} else
				map[i] = -1;
		}
		idx = 0;
		int new_adj[][] = new int[new_num_vet][3];
		for (int i = 0; i < v_num; i++) {
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

	public void toOrder(Params p, GeneOrder order) {

		for (int color = 0; color < 3; color++) {
			for (int i = 0; i < node_num; i++)
				unused[i] = true;

			int idx = 0;
			int gene;
			int current = -1, next;

			while (cap_vet[color].size() > 0) {
				order.gene_order[idx][color] = Const.CAP;
				idx++;
				next = cap_vet[color].pop_back();
				while (next != Const.CAP) {
					current = next;
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
					current = next;
					next = adj_mat[current][color];
				}
				cap_vet[color].remove(current);
				order.gene_order[idx][color] = Const.CAP;
				idx++;
				unused[current] = false;
			}

			for (short i = 0; i < node_num; i++) {
				if (!unused[i])
					continue;
				current = i;
				order.gene_order[idx][color] = Const.CAP;
				idx++;
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
				order.gene_order[idx][color] = Const.CAP;
				idx++;
			}
		}
		this.gene_order_to_graph(order);
	}
}
