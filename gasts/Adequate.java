package gasts;

// the class to detect adequate subgraphs on CMBG

// this does not support the option to avoid making circular chromosomes

import java.util.ArrayList;

public class Adequate {

	public static boolean AS1(GraphXu g, EF result) {
		boolean unused[] = new boolean[g.num_reg_vtx];
		for (int i = 0; i < g.num_reg_vtx; i++)
			unused[i] = true;
		for (int v = 0; v < g.num_reg_vtx; v++) {
			if (!unused[v])
				continue;
			int vn[] = g.neighbor[v];
			if (vn[0] < Constant.CAP1 && (vn[0] == vn[1] || vn[0] == vn[2])) {
				if (result.black.size() == 0)
					result.black.add(new ArrayList<Integer>());
				result.black.get(0).add(v);
				result.black.get(0).add(vn[0]);
				unused[vn[0]] = false;
			} else if (vn[1] == vn[2] && vn[1] < Constant.CAP1) {
				if (result.black.size() == 0)
					result.black.add(new ArrayList<Integer>());
				result.black.get(0).add(v);
				result.black.get(0).add(vn[1]);
				unused[vn[1]] = false;
				unused[vn[2]] = false;
			}
		}
		if (result.black.size() > 0) {
			result.info = "AS1";
			return true;
		} else
			return false;
	}

	public static boolean AS2(GraphXu g, EF result) {
		int four_cycle[] = new int[4];
		boolean two = false;
		boolean unused0[] = new boolean[g.num_reg_vtx];
		for (int i = 0; i < g.num_reg_vtx; i++)
			unused0[i] = true;

		for (int color = 0; color < 3; color++) {
			boolean unused[] = new boolean[g.num_reg_vtx];
			for (int i = 0; i < g.num_reg_vtx; i++)
				unused[i] = unused0[i];

			for (int i = 0; i < g.num_reg_vtx; i++) {
				if (!unused[i])
					continue;
				unused[i] = false;
				int j = g.neighbor[i][color];
				if (j >= Constant.CAP1 || !unused[j])
					continue;
				unused[j] = false;
				int c1 = (color + 1) % 3, c2 = (color + 2) % 3;
				int i1 = g.neighbor[i][c1], i2 = g.neighbor[i][c2];
				int j1 = g.neighbor[j][c1], j2 = g.neighbor[j][c2];
				if (i1 < Constant.CAP1 && j2 < Constant.CAP1
						&& g.neighbor[i1][color] == j2 && unused[i1]
						&& unused[j2]) {
					// System.out.println("it's e!");
					if (result.black.size() == 0)
						result.black.add(new ArrayList<Integer>());
					result.black.get(0).add(i);
					result.black.get(0).add(i1);
					result.black.get(0).add(j);
					result.black.get(0).add(j2);
					unused0[i] = unused0[j] = unused0[i1] = unused0[j2] = false;
					unused[i] = unused[j] = unused[i1] = unused[j2] = false;

					if (i == j || i == i1 || i == j2 || j == i1 || j == j2
							|| i1 == j2)
						System.out.println("error, two elements are the same!");
				} else if (i2 < Constant.CAP1 && j1 < Constant.CAP1
						&& g.neighbor[i2][color] == j1 && unused[i2]
						&& unused[j1]) {
					// System.out.println("it's f!");
					if (result.black.size() == 0)
						result.black.add(new ArrayList<Integer>());
					result.black.get(0).add(i);
					result.black.get(0).add(i2);
					result.black.get(0).add(j);
					result.black.get(0).add(j1);
					unused0[i] = unused0[j] = unused0[i2] = unused0[j1] = false;
					unused[i] = unused[j] = unused[i2] = unused[j1] = false;

					if (i == j || i == i2 || i == j1 || j == i2 || j == j1
							|| i2 == j1)
						System.out.println("error, two elements are the same!");
				} else if (i1 < Constant.CAP1
						&& i2 < Constant.CAP1
						&& // two four cycles
						j1 < Constant.CAP1 && j2 < Constant.CAP1
						&& g.neighbor[i1][color] == j1
						&& g.neighbor[i2][color] == j2 && i1 != j2 && i2 != j1
						&& unused[i1] && unused[i2] && unused[j1] && unused[j2]) {
					// System.out.println("it's g!");
					if (result.black.size() == 0)
						result.black.add(new ArrayList<Integer>());
					result.black.get(0).add(i);
					result.black.get(0).add(j);
					result.black.get(0).add(i1);
					result.black.get(0).add(j1);
					result.black.get(0).add(i2);
					result.black.get(0).add(j2);
					unused0[i] = unused0[j] = unused0[i1] = unused0[i2] = unused0[j1] = unused0[j2] = false;
					unused[i] = unused[j] = unused[i1] = unused[i2] = unused[j1] = unused[j2] = false;

					if (i == j || i == i1 || i == i2 || i == j1 || i == j2
							|| j == i1 || j == i2 || j == j1 || j == j2
							|| i1 == i2 || i1 == j1 || i1 == j2 || i2 == j1
							|| i2 == j2 || j1 == j2)
						System.out.println("error, two elements are the same!");
				} else if (i1 < Constant.CAP1 && j1 < Constant.CAP1
						&& g.neighbor[i1][color] == j1
						&& (i2 == j1 || j2 == i1) && unused[i1] && unused[j1]) {
					// System.out.println("it's a!");
					if (result.black.size() == 0)
						result.black.add(new ArrayList<Integer>());
					result.black.get(0).add(i);
					result.black.get(0).add(i1);
					result.black.get(0).add(j);
					result.black.get(0).add(j1);
					unused0[i] = unused0[j] = unused0[i1] = unused0[j1] = false;
					unused[i] = unused[j] = unused[i1] = unused[j1] = false;

					if (i == j || i == i1 || i == j1 || j == i1 || j == j1
							|| i1 == j1)
						System.out.println("error, two elements are the same!");
				} else if (i2 < Constant.CAP1 && j2 < Constant.CAP1
						&& g.neighbor[i2][color] == j2
						&& (i1 == j2 || j1 == i2) && unused[i2] && unused[j2]) {
					// System.out.println("it's b!");
					if (result.black.size() == 0)
						result.black.add(new ArrayList<Integer>());
					result.black.get(0).add(i);
					result.black.get(0).add(i2);
					result.black.get(0).add(j);
					result.black.get(0).add(j2);
					unused0[i] = unused0[j] = unused0[i2] = unused0[j2] = false;
					unused[i] = unused[j] = unused[i2] = unused[j2] = false;

					if (i == j || i == i2 || i == j2 || j == i2 || j == j2
							|| i2 == j2)
						System.out.println("error, two elements are the same!");
				} else if (result.black.size() == 0 && two == false) {
					if (i1 < Constant.CAP1 && j1 < Constant.CAP1
							&& g.neighbor[i1][color] == j1 && unused[i1]
							&& unused[j1]) {
						// System.out.println("it's c!");
						two = true;
						four_cycle[0] = i;
						four_cycle[1] = j;
						four_cycle[2] = i1;
						four_cycle[3] = j1;
						unused0[i] = unused0[j] = unused0[i1] = unused0[j1] = false;
						unused[i] = unused[j] = unused[i1] = unused[j1] = false;

						if (i == j || i == i1 || i == j1 || j == i1 || j == j1
								|| i1 == j1)
							System.out
									.println("error, two elements are the same!");
					} else if (i2 < Constant.CAP1 && j2 < Constant.CAP1
							&& g.neighbor[i2][color] == j2 && unused[i2]
							&& unused[j2]) {
						// System.out.println("it's d!");
						two = true;
						four_cycle[0] = i;
						four_cycle[1] = j;
						four_cycle[2] = i2;
						four_cycle[3] = j2;
						unused0[i] = unused0[j] = unused0[i2] = unused0[j2] = false;
						unused[i] = unused[j] = unused[i2] = unused[j2] = false;

						if (i == j || i == i2 || i == j2 || j == i2 || j == j2
								|| i2 == j2)
							System.out
									.println("error, two elements are the same!");
					}
				}
			}
		}
		if (result.black.size() > 0 || two) {
			if (result.black.size() > 0) {
				result.info = "AS2 m=1";
				return true;
			} else {
				result.black.add(new ArrayList<Integer>());
				result.black.add(new ArrayList<Integer>());
				result.black.get(0).add(four_cycle[0]);
				result.black.get(0).add(four_cycle[1]);
				result.black.get(0).add(four_cycle[2]);
				result.black.get(0).add(four_cycle[3]);
				result.black.get(1).add(four_cycle[0]);
				result.black.get(1).add(four_cycle[2]);
				result.black.get(1).add(four_cycle[1]);
				result.black.get(1).add(four_cycle[3]);
				result.info = "AS2 m=2";
				return true;
			}
		} else
			return false;
	}

	public static boolean AS4(GraphXu g, EF result) {
		result.black.clear();
		boolean unused[] = new boolean[g.num_reg_vtx];
		for (int i = 0; i < g.num_reg_vtx; i++)
			unused[i] = true;
		int temp[] = new int[8];
		for (int i = 0; i < g.num_reg_vtx; i++) {
			if (!unused[i])
				continue;
			outj: for (int j = (int) (i + 1); j < g.num_reg_vtx; j++) {
				if (!unused[j])
					continue;

				boolean flag = false;

				if (g.is_connected(i, j)) {
					// 5-3-5, 5-5-5, 3-5-5
					int color = g.incident(i, j);
					int c_u_l, c_u_r, c_d_l, c_d_r; // colors
					int p_u_l, p_u_r, p_d_l, p_d_r, p_u_m, p_d_m;
					c_u_l = c_d_r = (color + 1) % 3;
					c_d_l = c_u_r = (color + 2) % 3;
					p_u_l = g.neighbor[i][c_u_l];
					p_d_l = g.neighbor[i][c_d_l];
					p_u_r = g.neighbor[j][c_u_r];
					p_d_r = g.neighbor[j][c_d_r];
					if (p_u_l >= Constant.CAP1 || p_u_r >= Constant.CAP1
							|| p_d_l >= Constant.CAP1 || p_d_r >= Constant.CAP1
							|| !unused[p_u_l] || !unused[p_u_r]
							|| !unused[p_d_l] || !unused[p_d_r])
						continue;
					if (p_u_l == p_u_r) {
						// lower 3-5-5
						if (g.neighbor[p_d_l][c_d_r] != g.neighbor[p_d_r][c_d_l])
							continue;
						int m0 = g.neighbor[p_d_l][c_d_r];
						int m1 = g.neighbor[p_d_l][color];
						int m2 = g.neighbor[p_d_r][color];
						if (m0 >= Constant.CAP1 || m1 >= Constant.CAP1
								|| m2 >= Constant.CAP1 || !unused[m0]
								|| !unused[m1] || !unused[m2])
							continue;
						if (g.is_connected(m1, m2)) {
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
						if (g.neighbor[p_u_l][c_u_r] != g.neighbor[p_u_r][c_u_l])
							continue;
						int m0 = g.neighbor[p_u_l][c_u_r];
						int m1 = g.neighbor[p_u_l][color];
						int m2 = g.neighbor[p_u_r][color];
						if (m0 >= Constant.CAP1 || m1 >= Constant.CAP1
								|| m2 >= Constant.CAP1 || !unused[m0]
								|| !unused[m1] || !unused[m2])
							continue;
						if (g.is_connected(m1, m2)) {
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
							if (p_u_m >= Constant.CAP1
									|| p_d_m >= Constant.CAP1 || !unused[p_u_m]
									|| !unused[p_d_m])
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
					int left[] = new int[3];
					int right[] = new int[3];
					int count = 0;
					int remain = -1;
					for (int color = 0; color < 3; color++) {
						left[color] = g.neighbor[i][color];
						right[color] = g.neighbor[j][color];
						if (left[color] >= Constant.CAP1
								|| right[color] >= Constant.CAP1
								|| !unused[left[color]]
								|| !unused[right[color]])
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
					// System.out.println("3-3-3 or 3-3-6");
				} // end if-else i-j connected
				if (flag) {
					if (result.black.size() == 0)
						result.black.add(new ArrayList<Integer>());
					for (int k = 0; k < 8; k++) {
						result.black.get(0).add(temp[k]);
						unused[temp[k]] = false;
					}
				}
			} // end for j
		} // end for i
		if (result.black.size() > 0) {
			result.info += "AS4";
			return true;
		} else
			return false;
	} // end AS4

	public static void AS0(GraphXu g, EF result) {
		// assume g.num_reg_vtx>g.num_free_caps
		for (int i = 0; i < g.num_reg_vtx - 1; i++)
			result.black.add(new ArrayList<Integer>());
		for (int i = 1; i < g.num_reg_vtx; i++) {
			result.black.get(i - 1).add((Integer) (int) 0);
			result.black.get(i - 1).add((Integer) i);
		}
		if (g.num_free_caps > 0) {
			result.black.add(new ArrayList<Integer>());
			result.black.get(g.num_reg_vtx - 1).add((Integer) (int) 0);
			result.black.get(g.num_reg_vtx - 1).add((Integer) Constant.CAP0);
		}
		result.info = "AS0";
	}

	public static boolean C0AS1T(GraphXu g, EF result) {
		int num_free_caps = g.num_free_caps;

		for (int i = 0; i < g.num_reg_vtx && num_free_caps > 0; i++) {
			int count = 0;
			for (int color = 0; color < 3; color++) {
				if (g.neighbor[i][color] == Constant.CAP0)
					count++;
			}
			if (count >= 3) {
				num_free_caps--;
				if (result.black.size() == 0)
					result.black.add(new ArrayList<Integer>());
				result.black.get(0).add((Integer) i);
				result.black.get(0).add(Constant.CAP0);
			}
		}

		if (result.black.size() > 0) {
			result.info = "C0AS1-triple";
			return true;
		} else
			return false;
	}

	public static boolean C0AS1D(GraphXu g, EF result) {
		int num_free_caps = g.num_free_caps;
		for (int i = 0; i < g.num_reg_vtx && num_free_caps > 0; i++) {
			int count = 0;
			for (int color = 0; color < 3; color++) {
				if (g.neighbor[i][color] == Constant.CAP0)
					count++;
			}
			if (count >= 2) {
				if (result.black.size() == 0)
					result.black.add(new ArrayList<Integer>());
				result.black.get(0).add((Integer) i);
				result.black.get(0).add(Constant.CAP0);
				result.info = "C0AS1-double:ODD";
				return true;
			}
		}
		return false;
	}

	public static boolean C0AS23(GraphXu g, EF result) {
		Integer vertex[] = new Integer[3];
		for (int red = 0; red < g.cap0_adj.get(0).size(); red++) {
			vertex[0] = g.cap0_adj.get(0).get(red);
			for (int blue = 0; blue < g.cap0_adj.get(1).size(); blue++) {
				vertex[1] = g.cap0_adj.get(1).get(blue);
				if (vertex[1] == vertex[0]) {
					System.out.println("two vertices are identical");
					System.exit(-1);
					continue;
				}
				search: for (int green = 0; green < g.cap0_adj.get(2).size(); green++) {
					vertex[2] = g.cap0_adj.get(2).get(green);
					if (vertex[2] == vertex[0] || vertex[2] == vertex[1]) {
						System.out.println("two vertices are identical");
						System.exit(-1);
						continue;
					}

					int count = 0;
					int c_up = -1, c_down = -1;
					for (int i = 0; i < 2; i++) {
						for (int j = i + 1; j < 3; j++) {
							if (g.is_connected(vertex[i], vertex[j]))
								count++;
							else {
								c_up = i;
								c_down = j;
							}
						}
					}

					if (count == 3) {
						// find a size 2 CAS
						if (result.black.size() == 0)
							result.black.add(new ArrayList<Integer>());
						result.black.get(0).add(vertex[0]);
						result.black.get(0).add(Constant.CAP0);
						result.black.get(0).add(vertex[1]);
						result.black.get(0).add(vertex[2]);
						result.info = "C0AS23-2";
						return true;
					} else if (count == 2) {
						int color = 3 - c_up - c_down;
						int up = g.neighbor[vertex[c_up]][color];
						int down = g.neighbor[vertex[c_down]][color];
						if (g.is_connected(up, down)) {
							if (result.black.size() == 0)
								result.black.add(new ArrayList<Integer>());
							result.black.get(0).add(vertex[color]);
							result.black.get(0).add(Constant.CAP0);
							result.black.get(0).add(vertex[c_up]);
							result.black.get(0).add(vertex[c_down]);
							result.black.get(0).add(up);
							result.black.get(0).add(down);
							result.info = "C0AS23-3a";
							return true;
						}
					}

					for (int color = 0; color < 3; color++) {
						int up_mid = vertex[color];
						int c_left = (color + 1) % 3;
						int c_right = (color + 2) % 3;
						int down_left = vertex[c_left];
						int down_right = vertex[c_right];
						int up_left = g.neighbor[up_mid][c_left];
						int up_right = g.neighbor[up_mid][c_right];
						if (up_left == down_right || up_right == down_left)
							continue;
						if (up_left >= Constant.CAP1
								|| up_right >= Constant.CAP1)
							continue;
						if (g.is_connected(down_left, up_left)
								&& g.is_connected(down_right, up_right)) {
							if (result.black.size() == 0)
								result.black.add(new ArrayList<Integer>());
							result.black.get(0).add(up_mid);
							result.black.get(0).add(Constant.CAP0);
							result.black.get(0).add(up_left);
							result.black.get(0).add(down_left);
							result.black.get(0).add(up_right);
							result.black.get(0).add(down_right);
							result.info = "C0AS23-3b";
							// System.out.println("C0AS23-3b!!!");
							return true;
						}
					}

					for (int color0 = 0; color0 < 3; color0++) {
						for (int i = 1; i <= 2; i++) {
							int color1 = (color0 + i) % 3;
							int color2 = 3 - color0 - color1;

							if (g.neighbor[vertex[color1]][color0] == vertex[color2]) {
								int up_mid = g.neighbor[vertex[color0]][color1];
								int up_right = g.neighbor[vertex[color2]][color1];
								if (up_mid >= Constant.CAP1
										|| up_right >= Constant.CAP1)
									continue search;
								if (g.neighbor[vertex[0]][color2] == up_right
										&& g.neighbor[up_mid][color0] == up_right) {
									if (result.black.size() == 0)
										result.black
												.add(new ArrayList<Integer>());
									result.black.get(0).add(vertex[color0]);
									result.black.get(0).add(Constant.CAP0);
									result.black.get(0).add(vertex[color1]);
									result.black.get(0).add(up_mid);
									result.black.get(0).add(vertex[color2]);
									result.black.get(0).add(up_right);
									result.info = "C0AS23-3c";
									return true;
								}
							}
						}
					}
				}
			}
		}
		return false;
	}

	public static boolean C01AS1(GraphXu g, EF result) {
		boolean unused[] = new boolean[g.num_reg_vtx];
		for (int i = 0; i < g.num_reg_vtx; i++)
			unused[i] = true;

		for (int color = 0; color < 3; color++) {
			outc1: for (int c1 : g.cap1_adj.get(color)) {
				if (!unused[c1])
					continue;
				for (int c0 : g.cap0_adj.get(color)) {
					if (!unused[c0] || c1 == c0)
						continue;
					if (g.is_connected(c0, c1)) {
						if (result.black.size() == 0)
							result.black.add(new ArrayList<Integer>());
						result.black.get(0).add((Integer) c0);
						result.black.get(0).add((Integer) c1);
						unused[c0] = unused[c1] = false;
						continue outc1;
					}
				}
			}
		}
		if (result.black.size() == 0)
			return false;
		else {
			result.info = "C01AS1";
			// System.out.println("find a result in C01AS1");
			return true;
		}
	}

	public static boolean C01AS3(GraphXu g, EF result) {
		boolean unused[] = new boolean[g.num_reg_vtx];
		for (int i = 0; i < g.num_reg_vtx; i++)
			unused[i] = true;

		for (int color = 0; color < 3; color++) {
			int x1 = (color + 1) % 3;
			int x2 = (color + 2) % 3;
			outc1: for (int c1 : g.cap1_adj.get(color)) {
				if (!unused[c1])
					continue;
				int c1x1 = g.neighbor[c1][x1];
				int c1x2 = g.neighbor[c1][x2];
				if (c1x1 >= Constant.CAP1 || c1x2 >= Constant.CAP1
						|| !unused[c1x1] || !unused[c1x2])
					continue;
				for (int c0 : g.cap0_adj.get(color)) {
					if (!unused[c0] || c1 == c0)
						continue;
					int c0x1 = g.neighbor[c0][x1];
					int c0x2 = g.neighbor[c0][x2];
					if (c0x1 >= Constant.CAP1 || c0x2 >= Constant.CAP1
							|| !unused[c0x1] || !unused[c0x2])
						continue;
					if (g.is_connected(c0x1, c1x1)
							&& g.is_connected(c0x2, c1x2)) { // 3-3
						if (result.black.size() == 0)
							result.black.add(new ArrayList<Integer>());
						result.black.get(0).add((Integer) c0);
						result.black.get(0).add((Integer) c1);
						result.black.get(0).add((Integer) c0x1);
						result.black.get(0).add((Integer) c1x1);
						result.black.get(0).add((Integer) c0x2);
						result.black.get(0).add((Integer) c1x2);
						unused[c0] = unused[c1] = unused[c0x1] = unused[c0x2] = unused[c1x1] = unused[c1x2] = false;
						continue outc1;
					} // 3-3
					else if (g.is_connected(c0x1, c0x2)
							&& g.is_connected(c1x1, c1x2)
							&& (g.is_connected(c0x1, c1x1) || g.is_connected(
									c0x2, c1x2))) { // 3-3-1
						if (result.black.size() == 0)
							result.black.add(new ArrayList<Integer>());
						result.black.get(0).add((Integer) c0);
						result.black.get(0).add((Integer) c1);
						result.black.get(0).add((Integer) c0x1);
						result.black.get(0).add((Integer) c1x1);
						result.black.get(0).add((Integer) c0x2);
						result.black.get(0).add((Integer) c1x2);
						unused[c0] = unused[c1] = unused[c0x1] = unused[c0x2] = unused[c1x1] = unused[c1x2] = false;
						continue outc1;
					} // 3-3-1
					else if (c0x1 == c1x2 || c0x2 == c1x1) { // complex
						int triple1 = -1;
						int triple2 = -1;
						int double1 = -1;
						int double2 = -1;
						boolean flag = false;
						if (c0x1 == c1x2) {
							triple1 = c0x1;
							triple2 = g.neighbor[triple1][color];
							if (triple2 >= Constant.CAP1 || !unused[triple2])
								continue;
							if (triple2 == c0x2) {
								double2 = g.neighbor[triple2][x1];
								if (double2 >= Constant.CAP1
										|| !unused[double2])
									continue;
								if (g.is_connected(c1x1, double2)) {
									flag = true;
									double1 = c1x1;
								}
							} else if (triple2 == c1x1) {
								double2 = g.neighbor[triple2][x2];
								if (double2 >= Constant.CAP1
										|| !unused[double2])
									continue;
								if (g.is_connected(c0x2, double2)) {
									flag = true;
									double1 = c0x2;
								}
							}
						}

						else if (c0x2 == c1x1) {
							triple1 = c0x2;
							triple2 = g.neighbor[triple1][color];
							if (triple2 >= Constant.CAP1 || !unused[triple2])
								continue;
							if (triple2 == c0x1) {
								double2 = g.neighbor[triple2][x2];
								if (double2 >= Constant.CAP1
										|| !unused[double2])
									continue;
								if (g.is_connected(c1x2, double2)) {
									flag = true;
									double1 = c1x2;
								}
							} else if (triple2 == c1x2) {
								double2 = g.neighbor[triple2][x1];
								if (double2 >= Constant.CAP1
										|| !unused[double2])
									continue;
								if (g.is_connected(c0x1, double2)) {
									flag = true;
									double1 = c0x1;
								}
							}
						}

						if (flag) {
							if (result.black.size() == 0)
								result.black.add(new ArrayList<Integer>());
							result.black.get(0).add(c0);
							result.black.get(0).add(c1);
							result.black.get(0).add(triple1);
							result.black.get(0).add(triple2);
							result.black.get(0).add(double1);
							result.black.get(0).add(double2);
							unused[c0] = unused[c1] = unused[triple1] = unused[triple2] = unused[double1] = unused[double2] = false;
						}
					} // complex
				} // looping for c0
			} // looping for c1
		} // color
		if (result.black.size() == 0)
			return false;
		else {
			result.info = "C01AS3";
			// System.out.println("find a result in C01AS1");
			return true;
		}
	}

}
