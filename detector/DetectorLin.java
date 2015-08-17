package detector;

import graphs.Graph;
import tools.*;

public class DetectorLin extends Detector {

	@Override
	public void detect_ASs(Graph g, int up_num) {
		this.clean();
		if (AS1(g))
			return;
		else if (AS2_one(g)) {
			if (this.num_detected > 1) {
				this.transMajor();
				if (AS4(g))
					return;
				else if (C01AS1(g))
					return;
				else if (g.chr_num > 0 && C0AS1T(g))
					return;
				else if (g.chr_num > 0 && C0AS1D(g))
					return;
				else
					this.transMajorBack();
			}
		} else if (AS4(g))
			return;
		else if (C01AS1(g))
			return;
		else if (g.chr_num > 0 && C0AS1T(g))
			return;
		else if (g.chr_num > 0 && C0AS1D(g))
			return;
		else
			AS0(g);
		return;
	}

	public boolean AS1(Graph g) {
		for (int i = 0; i < g.v_num; i++) {
			if (g.check[i] == true)
				valid[i] = true;
			else
				valid[i] = false;
		}
		for (int v = 0; v < g.v_num; v++) {
			if (valid[v] == false)
				continue;
			vn[0] = g.adj_mat[v][0];
			vn[1] = g.adj_mat[v][1];
			vn[2] = g.adj_mat[v][2];
			if (vn[0] < Const.CUP && (vn[0] == vn[1] || vn[0] == vn[2])) {
				this.addVertex(v);
				this.addVertex(vn[0]);
				valid[vn[0]] = false;
			} else if (vn[1] == vn[2] && vn[1] < Const.CUP) {
				this.addVertex(v);
				this.addVertex(vn[1]);
				valid[vn[1]] = false;
				valid[vn[2]] = false;
			}
		}
		if (this.idx_major > 0) { // find AS1
			this.num_detected = 1;
			return true;
		} else
			return false;
	}

	// public boolean AS2_all(Graph g) {
	// // TODO Auto-generated method stub
	// return false;
	// }

	public boolean AS2_one(Graph g) {
		// TODO Auto-generated method stub
		boolean two = false;

		for (int i = 0; i < g.v_num; i++) {
			if (g.check[i] == true)
				valid0[i] = true;
			else
				valid0[i] = false;
		}

		for (int color = 0; color < 3; color++) {
			for (int i = 0; i < g.v_num; i++)
				valid[i] = valid0[i];

			for (int i = 0; i < g.v_num; i++) {
				if (valid[i] != true)
					continue;
				valid[i] = false;
				int j = g.adj_mat[i][color];
				if (j >= Const.CUP || valid[j] == false)
					continue;
				valid[j] = false;
				int c1 = (color + 1) % 3, c2 = (color + 2) % 3;
				int i1 = g.adj_mat[i][c1], i2 = g.adj_mat[i][c2];
				int j1 = g.adj_mat[j][c1], j2 = g.adj_mat[j][c2];
				if (i1 < Const.CUP && j2 < Const.CUP
						&& g.adj_mat[i1][color] == j2 && valid[i1] == true
						&& valid[j2] == true) {
					this.addVertex(i);
					this.addVertex(i1);
					this.addVertex(j);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i1] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i1] = valid[j2] = false;
				} else if (i2 < Const.CUP && j1 < Const.CUP
						&& g.adj_mat[i2][color] == j1 && valid[i2] == true
						&& valid[j1] == true) {
					this.addVertex(i);
					this.addVertex(i2);
					this.addVertex(j);
					this.addVertex(j1);
					valid0[i] = valid0[j] = valid0[i2] = valid0[j1] = false;
					valid[i] = valid[j] = valid[i2] = valid[j1] = false;
				} else if (i1 < Const.CUP
						&& i2 < Const.CUP
						&& // two four cycles
						j1 < Const.CUP && j2 < Const.CUP
						&& g.adj_mat[i1][color] == j1
						&& g.adj_mat[i2][color] == j2 && i1 != j2 && i2 != j1
						&& valid[i1] == true && valid[i2] == true
						&& valid[j1] == true && valid[j2] == true) {
					this.addVertex(i);
					this.addVertex(j);
					this.addVertex(i1);
					this.addVertex(j1);
					this.addVertex(i2);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i1] = valid0[i2] = valid0[j1] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i1] = valid[i2] = valid[j1] = valid[j2] = false;
				} else if (i1 < Const.CUP && j1 < Const.CUP
						&& g.adj_mat[i1][color] == j1 && (i2 == j1 || j2 == i1)
						&& valid[i1] == true && valid[j1] == true) {
					this.addVertex(i);
					this.addVertex(i1);
					this.addVertex(j);
					this.addVertex(j1);
					valid0[i] = valid0[j] = valid0[i1] = valid0[j1] = false;
					valid[i] = valid[j] = valid[i1] = valid[j1] = false;
				} else if (i2 < Const.CUP && j2 < Const.CUP
						&& g.adj_mat[i2][color] == j2 && (i1 == j2 || j1 == i2)
						&& valid[i2] == true && valid[j2] == true) {
					this.addVertex(i);
					this.addVertex(i2);
					this.addVertex(j);
					this.addVertex(j2);
					valid0[i] = valid0[j] = valid0[i2] = valid0[j2] = false;
					valid[i] = valid[j] = valid[i2] = valid[j2] = false;
				} else if (this.idx_major == 0 && two == false) {
					if (i1 < Const.CUP && j1 < Const.CUP
							&& g.adj_mat[i1][color] == j1 && valid[i1] == true
							&& valid[j1] == true) {
						two = true;
						four_cycle[0] = i;
						four_cycle[1] = j;
						four_cycle[2] = i1;
						four_cycle[3] = j1;
						valid0[i] = valid0[j] = valid0[i1] = valid0[j1] = false;
						valid[i] = valid[j] = valid[i1] = valid[j1] = false;
					} else if (i2 < Const.CUP && j2 < Const.CUP
							&& g.adj_mat[i2][color] == j2 && valid[i2] == true
							&& valid[j2] == true) {
						two = true;
						four_cycle[0] = i;
						four_cycle[1] = j;
						four_cycle[2] = i2;
						four_cycle[3] = j2;
						valid0[i] = valid0[j] = valid0[i2] = valid0[j2] = false;
						valid[i] = valid[j] = valid[i2] = valid[j2] = false;
					}
				}
			}
		}
		if (this.idx_major > 0 || two) {
			if (this.idx_major > 0) {
				this.num_detected = 1;
				return true;
			} else {
				this.addVertex(four_cycle[0]);
				this.addVertex(four_cycle[1]);
				this.addVertex(four_cycle[2]);
				this.addVertex(four_cycle[3]);
				this.addVertex(four_cycle[0]);
				this.addVertex(four_cycle[2]);
				this.addVertex(four_cycle[1]);
				this.addVertex(four_cycle[3]);
				this.num_detected = 2;
				return true;
			}
		} else
			return false;
	}

	public boolean AS4(Graph g) {
		// TODO Auto-generated method stub
		for (int i = 0; i < g.v_num; i++) {
			if (g.check[i] == true)
				valid[i] = true;
			else
				valid[i] = false;
		}

		for (int i = 0; i < g.v_num; i++) {
			if (valid[i] == false)
				continue;
			outj: for (int j = i + 1; j < g.v_num; j++) {
				// if (i == 94 && j == 106)
				// System.out.print("ok");
				if (valid[j] == false)
					continue;

				boolean flag = false;

				if (g.is_connected(i, j) == true) {
					// 5-3-5, 5-5-5, 3-5-5
					int color = g.incident(i, j);
					int c_u_l, c_u_r, c_d_l, c_d_r; // colors
					int p_u_l, p_u_r, p_d_l, p_d_r, p_u_m, p_d_m;
					c_u_l = c_d_r = (color + 1) % 3;
					c_d_l = c_u_r = (color + 2) % 3;
					p_u_l = g.adj_mat[i][c_u_l];
					p_d_l = g.adj_mat[i][c_d_l];
					p_u_r = g.adj_mat[j][c_u_r];
					p_d_r = g.adj_mat[j][c_d_r];
					if (p_u_l >= Const.CUP || p_u_r >= Const.CUP
							|| p_d_l >= Const.CUP || p_d_r >= Const.CUP
							|| valid[p_u_l] == false || valid[p_u_r] == false
							|| valid[p_d_l] == false || valid[p_d_r] == false)
						continue;
					if (p_u_l == p_u_r) {
						// lower 3-5-5
						if (g.adj_mat[p_d_l][c_d_r] != g.adj_mat[p_d_r][c_d_l])
							continue;
						int m0 = g.adj_mat[p_d_l][c_d_r];
						int m1 = g.adj_mat[p_d_l][color];
						int m2 = g.adj_mat[p_d_r][color];
						if (m0 >= Const.CUP || m1 >= Const.CUP
								|| m2 >= Const.CUP || valid[m0] == false
								|| valid[m1] == false || valid[m2] == false)
							continue;
						if (g.is_connected(m1, m2) == true) {
							// 3-5-5 l
							flag = true;
							temp[0] = i;
							temp[1] = p_d_l;
							temp[2] = j;
							temp[3] = p_d_r;
							temp[4] = p_u_l;
							temp[5] = m0;
							temp[6] = m1;
							temp[7] = m2;
							for (int ii = 0; ii < 8; ii++)
								for (int jj = ii + 1; jj < 8; jj++)
									if (temp[ii] == temp[jj])
										continue outj;
							// System.out.println("3-5-5l");
						} else
							continue;
					} else if (p_d_l == p_d_r) {
						// upper 3-5-5
						flag = true;
						if (g.adj_mat[p_u_l][c_u_r] != g.adj_mat[p_u_r][c_u_l])
							continue;
						int m0 = g.adj_mat[p_u_l][c_u_r];
						int m1 = g.adj_mat[p_u_l][color];
						int m2 = g.adj_mat[p_u_r][color];
						if (m0 >= Const.CUP || m1 >= Const.CUP
								|| m2 >= Const.CUP || valid[m0] == false
								|| valid[m1] == false || valid[m2] == false)
							continue;
						if (g.is_connected(m1, m2) == true) {
							// 3-5-5 l
							temp[0] = i;
							temp[1] = p_u_l;
							temp[2] = j;
							temp[3] = p_u_r;
							temp[4] = p_d_l;
							temp[5] = m0;
							temp[6] = m1;
							temp[7] = m2;
							for (int ii = 0; ii < 8; ii++)
								for (int jj = ii + 1; jj < 8; jj++)
									if (temp[ii] == temp[jj])
										continue outj;
							// System.out.println("3-5-5u");
						} else
							continue;
					} else {
						p_u_m = g.two_connected_by(p_u_l, p_u_r);
						p_d_m = g.two_connected_by(p_d_l, p_d_r);
						if (p_u_m != -1 && p_d_m != -1) {
							if (p_u_m >= Const.CUP || p_d_m >= Const.CUP
									|| valid[p_u_m] == false
									|| valid[p_d_m] == false)
								continue;
							// find two 5-cycles
							if (g.incident(p_u_l, p_d_l) != -1) {
								flag = true;
								temp[0] = i;
								temp[1] = p_d_l;
								temp[2] = j;
								temp[3] = p_u_l;
								temp[4] = p_d_r;
								temp[5] = p_d_m;
								temp[6] = p_u_r;
								temp[7] = p_u_m;
								for (int ii = 0; ii < 8; ii++)
									for (int jj = ii + 1; jj < 8; jj++)
										if (temp[ii] == temp[jj])
											continue outj;
								// System.out.println("5-3-5left");
							} else if (g.incident(p_u_r, p_d_r) != -1) {
								flag = true;
								temp[0] = i;
								temp[1] = p_d_r;
								temp[2] = j;
								temp[3] = p_u_r;
								temp[4] = p_d_l;
								temp[5] = p_d_m;
								temp[6] = p_u_l;
								temp[7] = p_u_m;
								for (int ii = 0; ii < 8; ii++)
									for (int jj = ii + 1; jj < 8; jj++)
										if (temp[ii] == temp[jj])
											continue outj;
								// System.out.println("5-3-5rght");
							} else if (g.incident(p_u_m, p_d_m) != -1) {
								flag = true;
								temp[0] = i;
								temp[1] = j;
								temp[2] = p_u_m;
								temp[3] = p_d_m;
								temp[4] = p_u_l;
								temp[5] = p_d_r;
								temp[6] = p_u_r;
								temp[7] = p_d_l;
								for (int ii = 0; ii < 8; ii++)
									for (int jj = ii + 1; jj < 8; jj++)
										if (temp[ii] == temp[jj])
											continue outj;
								// System.out.println("5-5-5");
							}
						} else
							continue;
					}
				} else {
					// 3-3-3, 3-3-6,

					int count = 0;
					int remain = -1;
					for (int color = 0; color < 3; color++) {
						left[color] = g.adj_mat[i][color];
						right[color] = g.adj_mat[j][color];
						if (left[color] >= Const.CUP
								|| right[color] >= Const.CUP
								|| valid[left[color]] == false
								|| valid[right[color]] == false)
							continue outj;
						if (g.incident(left[color], right[color]) != -1)
							count++;
						else
							remain = color;
					}

					if (count < 2)
						continue;
					else if (count == 2) {
						int l0 = left[remain];
						int r0 = right[remain];
						int l1 = left[(remain + 1) % 3];
						int r1 = right[(remain + 1) % 3];
						int l2 = left[(remain + 2) % 3];
						int r2 = right[(remain + 2) % 3];
						if (g.incident(l0, l1) != -1
								&& g.incident(r0, r1) != -1
								|| g.incident(l0, l2) != -1
								&& g.incident(r0, r2) != -1)
							;
						else
							continue;
					}
					// 3-3-3, 3-3-6
					flag = true;
					temp[0] = i;
					temp[1] = j;
					temp[2] = left[0];
					temp[3] = right[0];
					temp[4] = left[1];
					temp[5] = right[1];
					temp[6] = left[2];
					temp[7] = right[2];
					for (int ii = 0; ii < 8; ii++)
						for (int jj = ii + 1; jj < 8; jj++)
							if (temp[ii] == temp[jj])
								continue outj;
				} // end if-else i-j connected
				if (flag) {
					for (int k = 0; k < 8; k++) {
						this.addVertex(temp[k]);
						valid[temp[k]] = false;
					}
				}
			} // end for j
		} // end for i
		if (this.idx_major == 0)
			return false;
		this.num_detected = 1;
		return true;
	}

	public void AS0(Graph g) {
		int point = 0;
		if (!this.is_zero) {
			// use the heuristics to minimize the search space or maximize the
			// possible lower bound for heuristics
			g.get_start(0, 1);
			g.get_start(0, 2);
			g.get_start(1, 2);
			if (this.is_heu)
				point = g.cal_low_start();
			else
				point = g.cal_up_start();
		} else {
			// just select the point with lowest id
			for (int i = 0; i < g.v_num; i++)
				if (g.check[i]) {
					point = i;
					break;
				}
		}
		// add all possible combinations
		for (int i = 0; i < g.v_num; i++) {
			if (g.check[i] != false && i != point) {
				this.addVertex(point);
				this.addVertex(i);
			}
		}
		// if there are more chromosomes add more
		if (g.chr_num > 0) {
			this.addVertex(point);
			this.addVertex((int) Const.CAP);
		}
		this.num_detected = this.idx_major / 2;
	}

	public boolean C01AS1(Graph g) {
		for (int i = 0; i < g.v_num; i++) {
			if (g.check[i] == true)
				valid[i] = true;
			else
				valid[i] = false;
		}

		for (int color = 0; color < 3; color++) {
			outc1: for (int i = 0; i < g.cup_vet[color].size(); i++) {
				int c1 = g.cup_vet[color].c[i];
				if (valid[c1] == false)
					continue;
				for (int j = 0; j < g.cap_vet[color].size(); j++) {
					int c0 = g.cap_vet[color].c[j];
					if (valid[c0] == false || c1 == c0)
						continue;
					if (g.is_connected(c0, c1) == true) {
						this.addVertex(c0);
						this.addVertex(c1);
						valid[c0] = valid[c1] = false;
						continue outc1;
					}
				}
			}
		}
		if (this.idx_major == 0)
			return false;
		else {
			this.num_detected = 1;
			return true;
		}
	}

	public boolean C0AS1T(Graph g) {
		int c_size = g.chr_num;

		for (int i = 0; i < g.v_num && c_size > 0; i++) {
			if (g.check[i] == false)
				continue;
			int count = 0;
			for (int color = 0; color < 3; color++) {
				if (g.adj_mat[i][color] == Const.CAP)
					count++;
			}
			if (count >= 3) {
				c_size--;
				this.addVertex(i);
				this.addVertex(Const.CAP);
			}
		}

		if (this.idx_major > 0) {
			this.num_detected = 1;
			return true;
		} else
			return false;
	}

	public boolean C0AS1D(Graph g) {
		int c_size = g.chr_num;
		for (int i = 0; i < g.v_num && c_size > 0; i++) {
			if (g.check[i] == false)
				continue;
			int count = 0;
			for (int color = 0; color < 3; color++) {
				if (g.adj_mat[i][color] == Const.CAP)
					count++;
			}
			if (count >= 2) {
				this.addVertex(i);
				this.addVertex(Const.CAP);
				this.num_detected = 1;
				return true;
			}
		}
		return false;
	}

}
