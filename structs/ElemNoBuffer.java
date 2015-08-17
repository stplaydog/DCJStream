package structs;

import graphs.Graph;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import tools.*;
import scheduler.*;

import detector.Detector;

public class ElemNoBuffer extends Elem {
	public int[] parent;

	public int[] p_size;
	public int[] p_remain;

	public int[] child;

	public int[] c_size;

	BufferedWriter writer;

	public ElemNoBuffer(Info info, int ub, int thread_id) {
		super();
		// TODO Auto-generated constructor stub
		this.max_elem_num = info.max_elem_sz;
		this.is_expanded = false;
		this.elem_sz = max_elem_num * info.avg_node_num;
		this.elem_other_sz = max_elem_num + 1;
		this.ub = ub;
		info.th_total[thread_id][ub] = 0;
		this.is_first = true;
	}

	public void expand() {
		this.parent = new int[elem_sz];
		this.p_idx = 0;
		this.p_check = 0;
		this.p_size = new int[elem_other_sz];
		this.p_remain = new int[elem_other_sz];

		this.child = new int[elem_sz];
		this.c_size = new int[elem_other_sz];
		this.c_idx = 0;
		this.c_check = 0;
		this.is_expanded = true;
	}

	public void double_expand() {
		// expand parent
		int tmp_parent[] = new int[p_idx];
		for (int i = 0; i < p_idx; i++)
			tmp_parent[i] = parent[i];
		this.parent = new int[elem_sz * 2];
		for (int i = 0; i < p_idx; i++)
			parent[i] = tmp_parent[i];

		int tmp_p_size[] = new int[p_check];
		for (int i = 0; i < p_check; i++)
			tmp_p_size[i] = p_size[i];
		this.p_size = new int[elem_other_sz * 2];
		for (int i = 0; i < p_check; i++)
			p_size[i] = tmp_p_size[i];

		int tmp_p_remain[] = new int[p_check];
		for (int i = 0; i < p_check; i++)
			tmp_p_remain[i] = p_remain[i];
		this.p_remain = new int[elem_other_sz * 2];
		for (int i = 0; i < p_check; i++)
			p_remain[i] = tmp_p_remain[i];

		// expand child
		int tmp_child[] = new int[this.c_idx];
		for (int i = 0; i < this.c_idx; i++)
			tmp_child[i] = this.child[i];
		this.child = new int[elem_sz * 2];
		for (int i = 0; i < this.c_idx; i++)
			child[i] = tmp_child[i];

		int tmp_c_size[] = new int[this.c_check];
		for (int i = 0; i < this.c_check; i++)
			tmp_c_size[i] = this.c_size[i];
		this.c_size = new int[elem_sz * 2];
		for (int i = 0; i < this.c_check; i++)
			c_size[i] = tmp_c_size[i];

		// update status
		this.max_elem_num = max_elem_num * 2;
		this.elem_sz = this.elem_sz * 2;
		this.elem_other_sz = this.elem_other_sz * 2;
	}

	public void add_node(Graph g, Detector ade, int s, int e, Info info,
			int thread_id) {
		if (!this.is_expanded) {
			long start = System.currentTimeMillis();
			this.expand();
			long end = System.currentTimeMillis();
			info.vec_time[thread_id] += (end - start);
		}
		if (this.c_check >= this.max_elem_num) {
			if (this.ub == info.max_up[thread_id]) {
				if (this.max_elem_num == info.max_elem_sz) {
					long start = System.currentTimeMillis();
					this.double_expand();
					long end = System.currentTimeMillis();
					info.vec_time[thread_id] += (end - start);
				} else {
					System.out.println(c_check+" "+this.max_elem_num);
					this.half_toFile(thread_id, info);
				}
			} else
				this.toFile(thread_id, info);
		}
		// add parent
		long start = System.currentTimeMillis();
		if (this.is_first) {
			if ((this.p_idx + g.idx_tmp) >= this.elem_sz) {
				this.double_expand();
			}
			for (int i = 0; i < g.idx_tmp; i++) {
				this.parent[p_idx] = g.major_tmp[i];
				this.p_idx++;
			}
			this.is_first = false;
			this.p_size[p_check] = g.idx_tmp;
			this.p_remain[p_check] = 1;
			p_check++;
			info.actual_byte += g.idx_tmp * 4;
		} else
			this.p_remain[p_check - 1] += 1;
		// add child
		for (int i = s; i < e; i++) {
			this.child[c_idx] = ade.major[i];
			this.c_idx++;
		}
		this.c_size[c_check] = (e - s);
		c_check++;
		info.th_total[thread_id][ub]++;
		info.actual_byte += (e - s) * 4;
		long end = System.currentTimeMillis();
		info.vec_time[thread_id] += (end - start);
	}

	public boolean get_node(Graph g, Info info, int thread_id) {
		if (this.c_check == 0) {
			if (info.f_check[thread_id][this.ub] > 0)
				this.fromFile(thread_id, info);
			else if (LoadBalancer.balance(info, this.ub, thread_id))
				this.fromFile(thread_id, info);
			else
				return false;
		}
		// get parent
		long start = System.currentTimeMillis();
		g.idx_tmp = 0;
		int start_p = p_idx - p_size[p_check - 1];
		for (int i = start_p; i < p_idx; i++, g.idx_tmp++)
			g.major_tmp[g.idx_tmp] = parent[i];
		p_remain[p_check - 1]--;
		if (p_remain[p_check - 1] == 0) {
			p_check--;
			p_idx = start_p;
			info.actual_byte -= p_size[p_check] * 4;
		}
		// get child
		int start_c = c_idx - c_size[c_check - 1];
		for (int i = start_c; i < c_idx; i++, g.idx_tmp++)
			g.major_tmp[g.idx_tmp] = child[i];
		c_idx = start_c;
		c_check--;
		info.actual_byte -= c_size[c_check] * 4;
		info.th_total[thread_id][ub]--;
		long end = System.currentTimeMillis();
		info.vec_time[thread_id] += (end - start);
		return true;
	}

	public void fork(int gran, Elem elem, int from_thread, int to_thread,
			Info info) {
		long start = System.currentTimeMillis();
		this.add_fork_node(gran, (ElemNoBuffer) elem, from_thread, to_thread,
				info);
		long end = System.currentTimeMillis();
		info.vec_time[from_thread] += (end - start);
	}

	public void add_fork_node(int gran, ElemNoBuffer col, int from_thread,
			int to_thread, Info info) {
		if (!this.is_expanded)
			this.expand();
		int p[] = new int[this.max_elem_num];
		int c[] = new int[this.max_elem_num];
		int p_len = 0;
		int c_len = 0;

		for (int num = 0; num < gran; num++) {
			if (num > this.max_elem_num)
				this.expand();
			// get parent
			p_len = 0;
			int start_p = col.p_idx - col.p_size[col.p_check - 1];
			for (int i = start_p; i < col.p_idx; i++, p_len++)
				p[p_len] = col.parent[i];
			col.p_remain[col.p_check - 1]--;
			// get child
			c_len = 0;
			int start_c = col.c_idx - col.c_size[col.c_check - 1];
			for (int i = start_c; i < col.c_idx; i++, c_len++)
				c[c_len] = col.child[i];
			col.c_idx = start_c;
			col.c_check--;
			info.th_total[from_thread][col.ub]--;
			// /////////////////////////
			// add parent
			if (this.is_first) {
				if ((this.p_idx + p_len) >= this.elem_sz) {
					this.double_expand();
				}
				for (int i = 0; i < p_len; i++) {
					this.parent[p_idx] = p[i];
					this.p_idx++;
				}
				this.is_first = false;
				this.p_size[p_check] = p_len;
				this.p_remain[p_check] = 1;
				p_check++;
			} else {
				if (this.is_first == false && p_check == 0)
					return;
				this.p_remain[p_check - 1] += 1;
			}
			// add child
			if ((this.c_idx + c_len) >= this.elem_sz) {
				this.double_expand();
			}
			for (int i = 0; i < c_len; i++) {
				this.child[c_idx] = c[i];
				this.c_idx++;
			}
			this.c_size[c_check] = c_len;
			c_check++;
			info.th_total[to_thread][ub]++;
			if (col.p_remain[col.p_check - 1] == 0) {
				col.p_check--;
				col.p_idx = start_p;
				this.is_first = true;
			}
		}
		this.refresh();
	}

	public void toFile(int thread_id, Info info) {
		long start_io = System.currentTimeMillis();
		for (int i = 0; i < 100; i++) {
			if (info.is_lock[thread_id])
				break;
			for (int j = 0; j < 100000; j++)
				;
		}
		info.f_check[thread_id][ub] += 1;

		String name = info.root_dir + "/" + thread_id + "_" + ub + "_"
				+ info.f_check[thread_id][ub];
		File file = new File(name);
		try {
			writer = new BufferedWriter(new FileWriter(file));

			for (int i = 0; i < p_idx; i++) {
				writer.write(parent[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = 0; i < p_check; i++) {
				writer.write(p_remain[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = 0; i < p_check; i++) {
				writer.write(p_size[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = 0; i < c_idx; i++) {
				writer.write(child[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = 0; i < c_check; i++) {
				writer.write(c_size[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();
			writer.flush();

			this.p_idx = 0;
			this.p_check = 0;
			this.c_idx = 0;
			this.c_check = 0;
			writer.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long end_io = System.currentTimeMillis();
		info.io_time[thread_id] += (end_io - start_io);
		info.space_usage[thread_id][Const.DSC_USG] += (float) file.length()
				/ Const.MB;
		this.refresh();
	}

	public void half_toFile(int thread_id, Info info) {

		long start_io = System.currentTimeMillis();
		for (int i = 0; i < 100; i++) {
			if (info.is_lock[thread_id])
				break;
			for (int j = 0; j < 100000; j++)
				;
		}
		info.f_check[thread_id][ub] += 1;

		String name = info.root_dir + "/" + thread_id + "_" + ub + "_"
				+ info.f_check[thread_id][ub];
		File file = new File(name);

		int total = 0;
		int tmp_p_idx = p_idx;
		int tmp_p_check = p_check;
		while (total < this.max_elem_num) {
			tmp_p_idx -= p_size[tmp_p_check - 1];
			total += p_remain[tmp_p_check - 1];
			tmp_p_check--;
		}

		int c_total = 0;
		int tmp_c_idx = c_idx;
		int tmp_c_check = c_check;
		while (c_total < total) {
			tmp_c_idx -= c_size[tmp_c_check - 1];
			c_total++;
			tmp_c_check--;
		}

		try {
			writer = new BufferedWriter(new FileWriter(file));

			for (int i = tmp_p_idx; i < p_idx; i++) {
				writer.write(parent[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = tmp_p_check; i < p_check; i++) {
				writer.write(p_remain[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = tmp_p_check; i < p_check; i++) {
				writer.write(p_size[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = tmp_c_idx; i < c_idx; i++) {
				writer.write(child[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();

			for (int i = tmp_c_check; i < c_check; i++) {
				writer.write(c_size[i] + " ");
				if (i % 10000 == 0 && i > 0)
					writer.newLine();
			}
			if ((p_idx - 1) % 10000 != 0)
				writer.newLine();
			writer.write("next");
			writer.newLine();
			writer.flush();

			this.p_idx = tmp_p_idx;
			this.p_check = tmp_p_check;
			this.c_idx = tmp_c_idx;
			this.c_check = tmp_c_check;
			writer.close();
			System.exit(1);
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		long end_io = System.currentTimeMillis();
		info.io_time[thread_id] += (end_io - start_io);
		info.space_usage[thread_id][Const.DSC_USG] += (float) file.length()
				/ Const.MB;
		this.refresh();

	}

	public boolean fromFile(int thread_id, Info info) {
		long start_io = System.currentTimeMillis();
		for (int i = 0; i < 100; i++) {
			if (info.is_lock[thread_id])
				break;
			for (int j = 0; j < 100000; j++)
				;
		}
		if (info.f_check[thread_id][ub] == 0)
			return false;

		String name = info.root_dir + "/" + thread_id + "_" + ub + "_"
				+ info.f_check[thread_id][ub];
		File file = new File(name);
		info.space_usage[thread_id][Const.DSC_USG] -= (float) file.length()
				/ Const.MB;

		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line[] = reader.readLine().trim().split(" ");
			while (line[0].indexOf("next") == -1) {
				if ((this.p_idx + line.length) > this.elem_sz)
					this.double_expand();
				for (int i = 0; i < line.length; i++) {
					parent[this.p_idx] = Integer.parseInt(line[i]);
					this.p_idx++;
				}
				line = reader.readLine().trim().split(" ");
			}
			line = reader.readLine().trim().split(" ");
			while (line[0].indexOf("next") == -1) {
				for (int i = 0; i < line.length; i++) {
					p_remain[this.p_check] = Integer.parseInt(line[i]);
					this.p_check++;
				}
				line = reader.readLine().trim().split(" ");
			}
			line = reader.readLine().trim().split(" ");
			this.p_check = 0;
			while (line[0].indexOf("next") == -1) {
				for (int i = 0; i < line.length; i++) {
					p_size[this.p_check] = Integer.parseInt(line[i]);
					this.p_check++;
				}
				line = reader.readLine().trim().split(" ");
			}

			line = reader.readLine().trim().split(" ");
			while (line[0].indexOf("next") == -1) {
				if (this.c_idx + line.length > this.elem_sz)
					this.double_expand();
				for (int i = 0; i < line.length; i++) {
					child[this.c_idx] = Integer.parseInt(line[i]);
					this.c_idx++;
				}
				line = reader.readLine().trim().split(" ");
			}

			line = reader.readLine().trim().split(" ");
			while (line[0].indexOf("next") == -1) {
				for (int i = 0; i < line.length; i++) {
					c_size[this.c_check] = Integer.parseInt(line[i]);
					this.c_check++;
				}
				line = reader.readLine().trim().split(" ");
			}

			// System.out.println("p_idx: " + this.p_idx + " p_check: "
			// + this.p_check + " c_idx: " + this.c_idx + " c_check: "
			// + this.c_check);
			reader.close();
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		info.f_check[thread_id][ub]--;
		file.delete();
		long end_io = System.currentTimeMillis();
		info.io_time[thread_id] += (end_io - start_io);

		return true;
	}

	public void p_r()
	{
		int start_p = p_idx - p_size[p_check - 1];
		for (int i = start_p; i < p_idx; i++)
			System.out.print(parent[i]+" ");
		System.out.println(p_remain[p_check-1]);
		System.out.println(p_size[p_check-1]);
	}

}
