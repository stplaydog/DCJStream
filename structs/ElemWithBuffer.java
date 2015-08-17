package structs;

import graphs.Graph;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import tools.Func;
import tools.Info;
import detector.Detector;

public class ElemWithBuffer extends Elem {
	public char[] write_buf;

	public int[][] parent;

	public int[][] p_size;
	public int[][] p_remain;

	public int[][] child;

	public int[][] c_size;

	boolean is_avail[];
	boolean is_full[];

	int final_p_idx[];
	int final_p_check[];
	int final_c_idx[];
	int final_c_check[];
	int time_stamp[];

	int cur_pos;
	int buff_sz;

	int global_time;

	FileWriter writer;

	public ElemWithBuffer(Info info, int ub, int thread_id) {
		super();
		// TODO Auto-generated constructor stub
		this.max_elem_num = info.max_elem_sz;
		this.is_expanded = false;
		this.elem_sz = max_elem_num * 20;
		this.elem_other_sz = max_elem_num + 1;
		this.ub = ub;
		info.th_total[thread_id][ub] = 0;
		this.is_first = true;
		this.cur_pos = 0;
		this.buff_sz = 4;
	}

	public void expand() {
		this.parent = new int[buff_sz][elem_sz];
		this.p_idx = 0;
		this.p_check = 0;
		this.p_size = new int[buff_sz][elem_other_sz];
		this.p_remain = new int[buff_sz][elem_other_sz];
		this.write_buf = new char[elem_sz * 2];

		this.child = new int[buff_sz][elem_sz];
		this.c_size = new int[buff_sz][elem_other_sz];

		this.is_avail = new boolean[buff_sz];
		this.is_full = new boolean[buff_sz];

		final_p_idx = new int[buff_sz];
		final_p_check = new int[buff_sz];
		final_c_idx = new int[buff_sz];
		final_c_check = new int[buff_sz];
		this.time_stamp = new int[buff_sz];

		for (int i = 0; i < buff_sz; i++) {
			is_avail[i] = true;
			is_full[i] = false;
		}

		this.c_idx = 0;
		this.c_check = 0;
		this.is_expanded = true;
		this.global_time = 0;

		for (int i = 0; i < this.buff_sz; i++)
			this.time_stamp[i] = 0;
	}

	public void add_node(Graph g, Detector ade, int s, int e, Info info,
			int thread_id) {
		if (!this.is_expanded)
			this.expand();

		if (this.c_check >= info.max_elem_sz) {
			this.switch_write_pos(info, thread_id);
		}

		// add parent
		if (this.is_first) {
			for (int i = 0; i < g.idx_tmp; i++) {
				this.parent[this.cur_pos][p_idx] = g.major_tmp[i];
				this.p_idx++;
			}
			this.is_first = false;
			this.p_size[this.cur_pos][p_check] = g.idx_tmp;
			this.p_remain[this.cur_pos][p_check] = 1;
			p_check++;
		} else
			this.p_remain[this.cur_pos][p_check - 1] += 1;
		// add child
		for (int i = s; i < e; i++) {
			this.child[this.cur_pos][c_idx] = ade.major[i];
			this.c_idx++;
		}
		this.c_size[this.cur_pos][c_check] = (e - s);
		c_check++;
		info.th_total[thread_id][ub]++;
	}

	public void switch_write_pos(Info info, int thread_id) {
		to_idx();
		this.cur_pos = detect_avail_write_pos();
		this.time_stamp[this.cur_pos] = (++this.global_time);
		if (this.amount_full() >= 2) {
			long start_th = System.currentTimeMillis();
			FileThread f_t = new FileThread(thread_id, info, this, true);
			long end_th = System.currentTimeMillis();
			info.thread_time[thread_id] += (end_th - start_th);
			f_t.run();
		}
		this.refresh();
	}

	public void to_idx() {
		is_full[cur_pos] = true;
		final_p_idx[cur_pos] = p_idx;
		final_p_check[cur_pos] = p_check;
		final_c_idx[cur_pos] = c_idx;
		final_c_check[cur_pos] = c_check;

		p_idx = 0;
		p_check = 0;
		c_idx = 0;
		c_check = 0;
	}

	public int detect_avail_write_pos() {
		int time = 10000000;
		int pos = -1;
		for (int i = 0; i < this.buff_sz; i++)
			if (!(this.is_full[i]) && this.is_avail[i])
				if (this.time_stamp[i] < time) {
					time = this.time_stamp[i];
					pos = i;
				}
		return pos;
	}

	public boolean get_node(Graph g, Info info, int thread_id) {
		if (this.c_check == 0)
			if (!this.switch_read_pos(info, thread_id))
				return false;
		// get parent
		g.idx_tmp = 0;
		int start_p = p_idx - p_size[this.cur_pos][p_check - 1];
		for (int i = start_p; i < p_idx; i++, g.idx_tmp++)
			g.major_tmp[g.idx_tmp] = parent[this.cur_pos][i];
		p_remain[this.cur_pos][p_check - 1]--;
		if (p_remain[this.cur_pos][p_check - 1] == 0) {
			p_check--;
			p_idx = start_p;
		}
		// get child
		int start_c = c_idx - c_size[this.cur_pos][c_check - 1];
		for (int i = start_c; i < c_idx; i++, g.idx_tmp++)
			g.major_tmp[g.idx_tmp] = child[this.cur_pos][i];
		c_idx = start_c;
		c_check--;
		info.th_total[thread_id][ub]--;
		return true;
	}

	public boolean switch_read_pos(Info info, int thread_id) {
		this.cur_pos = detect_avail_read_pos();
		if (this.cur_pos != -1) {
			is_full[cur_pos] = false;
			p_idx = final_p_idx[cur_pos];
			p_check = final_p_check[cur_pos];
			c_idx = final_c_idx[cur_pos];
			c_check = final_c_check[cur_pos];
		} else {
			if (info.f_check[thread_id][ub] == 0)
				return false;
			else {
				this.cur_pos = 0;
				// at this step, has to be synchronized.
				this.fromFile(thread_id, info);
			}
		}

		if (this.amount_full() < 2) {
			// read files
			long start_th = System.currentTimeMillis();
			FileThread f_t = new FileThread(thread_id, info, this, false);
			f_t.run();
			long end_th = System.currentTimeMillis();
			info.thread_time[thread_id] += (end_th - start_th);
		}
		return true;
	}

	public int detect_avail_read_pos() {
		int time = 0;
		int pos = -1;
		for (int i = 0; i < this.buff_sz; i++)
			if ((this.is_full[i]) && this.is_avail[i])
				if (this.time_stamp[i] > time) {
					time = this.time_stamp[i];
					pos = i;
				}
		return pos;
	}

	public int amount_full() {
		int result = 0;
		for (int i = 0; i < this.buff_sz; i++)
			if (this.is_full[i])
				result++;
		return result;
	}

	public void fork(int gran, Elem elem, int from_thread, int to_thread,
			Info info) {
		this.add_fork_node(gran, (ElemWithBuffer) elem, from_thread, to_thread,
				info);
	}

	public void add_fork_node(int gran, ElemWithBuffer col, int from_thread,
			int to_thread, Info info) {
		if (!this.is_expanded)
			this.expand();
		int p[] = new int[this.max_elem_num];
		int c[] = new int[this.max_elem_num];
		int p_len = 0;
		int c_len = 0;

		for (int num = 0; num < gran; num++) {
			// get parent
			p_len = 0;
			int start_p = col.p_idx - col.p_size[col.cur_pos][col.p_check - 1];
			for (int i = start_p; i < col.p_idx; i++, p_len++)
				p[p_len] = col.parent[col.cur_pos][i];
			col.p_remain[col.cur_pos][col.p_check - 1]--;
			// get child
			c_len = 0;
			int start_c = col.c_idx - col.c_size[col.cur_pos][col.c_check - 1];
			for (int i = start_c; i < col.c_idx; i++, c_len++)
				c[c_len] = col.child[col.cur_pos][i];
			col.c_idx = start_c;
			col.c_check--;
			info.th_total[from_thread][col.ub]--;
			// /////////////////////////
			// add parent
			if (this.is_first) {
				for (int i = 0; i < p_len; i++) {
					this.parent[this.cur_pos][p_idx] = p[i];
					this.p_idx++;
				}
				this.is_first = false;
				this.p_size[this.cur_pos][p_check] = p_len;
				this.p_remain[this.cur_pos][p_check] = 1;
				p_check++;
			} else
				this.p_remain[this.cur_pos][p_check - 1] += 1;
			// add child
			for (int i = 0; i < c_len; i++) {
				this.child[this.cur_pos][c_idx] = c[i];
				this.c_idx++;
			}
			this.c_size[this.cur_pos][c_check] = c_len;
			c_check++;
			info.th_total[to_thread][ub]++;
			if (col.p_remain[col.cur_pos][col.p_check - 1] == 0) {
				col.p_check--;
				col.p_idx = start_p;
				this.is_first = true;
			}
		}
		this.refresh();
	}

	public void toFile(int thread_id, Info info) {
		while (info.is_lock[thread_id])
			for (int i = 0; i < 10000; i++)
				;

		int w_pos = 0;
		int time = 100000;

		for (int i = 0; i < this.buff_sz; i++)
			if (this.is_full[i] && i != this.cur_pos
					&& this.time_stamp[i] < time) {
				w_pos = i;
				time = this.time_stamp[i];
			}
		this.is_avail[w_pos] = false;
		info.f_check[thread_id][ub] += 1;

		String name = info.root_dir + "/" + thread_id + "_" + ub + "_"
				+ info.f_check[thread_id][ub];

		File file = new File(name);

		try {
			long start_io = System.currentTimeMillis();
			FileWriter writer = new FileWriter(file);

			for (int i = 0; i < final_p_idx[w_pos]; i++)
				writer.write(parent[w_pos][i] + " ");
			writer.write("\n");

			// int len = this.to_write_buf(parent[w_pos], final_p_idx[w_pos]);
			// writer.write(this.write_buf, 0, len);
			// writer.write("\n");

			for (int i = 0; i < final_p_check[w_pos]; i++)
				writer.write(p_remain[w_pos][i] + " ");
			writer.write("\n");
			// len = this.to_write_buf(p_remain[w_pos], final_p_check[w_pos]);
			// writer.write(this.write_buf, 0, len);
			// writer.write("\n");
			for (int i = 0; i < final_p_check[w_pos]; i++)
				writer.write(p_size[w_pos][i] + " ");
			writer.write("\n");
			// len = this.to_write_buf(p_size[w_pos], final_p_check[w_pos]);
			// writer.write(this.write_buf, 0, len);
			// writer.write("\n");

			for (int i = 0; i < final_c_idx[w_pos]; i++)
				writer.write(child[w_pos][i] + " ");
			writer.write("\n");

			// len = this.to_write_buf(child[w_pos], final_c_idx[w_pos]);
			// writer.write(this.write_buf, 0, len);
			// writer.write("\n");
			for (int i = 0; i < final_c_check[w_pos]; i++)
				writer.write(c_size[w_pos][i] + " ");
			writer.write("\n");
			// len = this.to_write_buf(c_size[w_pos], final_c_check[w_pos]);
			// writer.write(this.write_buf, 0, len);
			// writer.write("\n");

			// writer.flush();
			writer.close();
			long end_io = System.currentTimeMillis();
			info.io_time[thread_id] += (end_io - start_io);

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		this.refresh();

		this.final_p_idx[w_pos] = 0;
		this.final_p_check[w_pos] = 0;
		this.final_c_idx[w_pos] = 0;
		this.final_c_check[w_pos] = 0;
		this.is_avail[w_pos] = true;
		this.is_full[w_pos] = false;
	}

	public void p_r(){}
	public int to_write_buf(int content[], int len) {
		int j = 0;
		for (int i = 0; i < len; i++) {
			char res[] = Func.itoa(content[i]);
			for (int k = 0; k < res.length; k++) {
				if (res[k] == '\0')
					break;
				this.write_buf[j] = res[k];
				j++;
			}
			this.write_buf[j] = ' ';
			j++;
		}
		return (j - 1);
	}

	public boolean fromFile(int thread_id, Info info) {
		while (info.is_lock[thread_id])
			for (int i = 0; i < 10000; i++)
				;

		if (info.f_check[thread_id][ub] == 0)
			return false;

		int r_pos = 0;
		int time = 100000;
		for (int i = 0; i < this.buff_sz; i++)
			if (!this.is_full[i] && i != this.cur_pos
					&& this.time_stamp[i] < time) {
				r_pos = i;
				time = this.time_stamp[i];
			}
		this.is_avail[r_pos] = false;
		// long start_io = System.currentTimeMillis();
		String name = info.root_dir + "/" + thread_id + "_" + ub + "_"
				+ info.f_check[thread_id][ub];
		File file = new File(name);
		try {
			BufferedReader reader = new BufferedReader(new FileReader(file));
			String line[] = reader.readLine().trim().split(" ");
			this.final_p_idx[r_pos] = line.length;
			for (int i = 0; i < final_p_idx[r_pos]; i++)
				parent[r_pos][i] = Integer.parseInt(line[i]);
			line = reader.readLine().trim().split(" ");
			this.final_p_check[r_pos] = line.length;
			for (int i = 0; i < final_p_check[r_pos]; i++)
				p_remain[r_pos][i] = Integer.parseInt(line[i]);
			line = reader.readLine().trim().split(" ");
			for (int i = 0; i < final_p_check[r_pos]; i++)
				p_size[r_pos][i] = Integer.parseInt(line[i]);

			line = reader.readLine().trim().split(" ");
			this.final_c_idx[r_pos] = line.length;
			for (int i = 0; i < final_c_idx[r_pos]; i++)
				child[r_pos][i] = Integer.parseInt(line[i]);
			line = reader.readLine().trim().split(" ");
			this.final_c_check[r_pos] = line.length;
			for (int i = 0; i < final_c_check[r_pos]; i++)
				c_size[r_pos][i] = Integer.parseInt(line[i]);

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
		// long end_io = System.currentTimeMillis();
		// info.io_time[thread_id] += (end_io - start_io);

		this.is_avail[r_pos] = true;
		this.is_full[r_pos] = true;
		return true;
	}

	class FileThread implements Runnable {
		int thread_id;
		Info info;
		ElemWithBuffer elem;
		boolean is_write;

		public FileThread(int thread_id, Info info, ElemWithBuffer elem,
				boolean is_write) {
			super();
			this.thread_id = thread_id;
			this.info = info;
			this.elem = elem;
			this.is_write = is_write;
		}

		@Override
		public void run() {
			// TODO Auto-generated method stub
			if (!this.is_write)
				elem.fromFile(thread_id, info);
			else
				elem.toFile(thread_id, info);
		}

	}

}
