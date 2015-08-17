package graphs;

import java.util.ArrayList;
import java.util.Random;

public class CList {
	public int c[];
	public int idx;
	int max_num;

	public CList() {
		super();
		this.idx = 0;
		this.max_num = 0;
	}

	public void init(int max_num) {
		this.max_num = max_num;
		this.c = new int[max_num];
		this.idx = -1;
	}

	public void init(CList l) {
		this.idx = l.idx;
		for (int i = 0; i <= l.idx; i++)
			this.c[i] = l.c[i];
	}

	public void init(ArrayList<Integer> cap_adj) {
		this.idx = cap_adj.size()-1;
		for (int i = 0; i < cap_adj.size(); i++)
			this.c[i] = cap_adj.get(i);

	}

	public boolean is_in(int query) {
		int start = 0;
		int end = idx;
		while (end >= start) {
			int current = (end - start) / 2;
			if (c[current] == query)
				return true;
			else if (current > query)
				end = current;
			else
				start = current;
		}
		return false;
	}

	public int pop_back() {
		int result = c[idx];
		c[idx] = 0;
		idx--;
		return result;
	}

	public void add(int elem) {
		idx++;
		if (idx == 0) {
			c[idx] = elem;
			return;
		}
		for (int i = 0; i < idx; i++) {
			if (i == 0 && c[i] > elem) {
				move_back(0);
				c[0] = elem;
				return;
			} else if (i != 0 && c[i - 1] < elem && c[i] > elem) {
				move_back(i);
				c[i] = elem;
				return;
			}
		}
		c[idx] = elem;
	}

	public void move_back(int pos) {
		for (int i = idx - 1; i >= pos; i--)
			c[i + 1] = c[i];
	}

	public void remove(int query) {
		int start = 0;
		int end = idx;
		while (end >= start) {
			int current = start
					+ (int) Math.ceil(((double) end - (double) start) / 2.0);
			if (c[current] == query) {
				move_forward(current);
				idx--;
				return;
			} else if (c[current] > query)
				end = current - 1;
			else
				start = current + 1;
		}
	}

	public void move_forward(int pos) {
		for (int i = pos; i < idx; i++)
			c[i] = c[i + 1];
		c[idx] = 0;
	}

	public void clear() {
		idx = -1;
	}

	public int size() {
		return idx + 1;
	}

	public static void main(String args[]) {
		Random rnd = new Random();
		CList list = new CList();
		list.init(1000);
		for (int i = 0; i < 100; i++) {
			int elem = i + 1;
			list.add(elem);
		}

		for (int i = 1; i < 10; i++) {
			list.remove(list.c[i * 10]);
		}

		for (int i = 0; i < 10; i++) {
			System.out.println(list.pop_back());
		}
	}
}
