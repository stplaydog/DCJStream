package detector;

import graphs.Graph;

public abstract class Detector {
	public int major[];
	int major_tmp[];
	boolean valid[];
	boolean valid0[];
	int idx_major;
	int idx_tmp;
	int num_detected = 0;
	int num_detected_tmp = 0;
	int four_cycle[];
	int[] tria;
	int start[][];
	int[] pI; // vertices of distance 1
	int[][] pII;// vertices of distance 2
	int[][] pIII;// vertices of distance 3
	int[] po;
	int vn[];
	int temp[];
	int left[];
	int right[];
	boolean is_heu;
	boolean is_zero;

	public long as_time;

	public abstract void detect_ASs(Graph g, int up_num);

	public void addVertex(int v) {
		major[idx_major] = v;
		idx_major++;
	}

	public void transMajor() {
		for (int i = 0; i < idx_major; i++)
			major_tmp[i] = major[i];
		idx_tmp = idx_major;
		idx_major = 0;
		num_detected_tmp = num_detected;
		num_detected = 0;
	}

	public void transMajorBack() {
		for (int i = 0; i < idx_tmp; i++)
			major[i] = major_tmp[i];
		idx_major = idx_tmp;
		num_detected = num_detected_tmp;
		num_detected_tmp = 0;
		idx_tmp = 0;
	}

	public int getNum_detected() {
		return num_detected;
	}

	public void setNum_detected(int num_detected) {
		this.num_detected = num_detected;
	}

	public int getIdx_major() {
		return idx_major;
	}

	public void setIdx_major(int idx_major) {
		this.idx_major = idx_major;
	}

	public void clean() {
		idx_major = 0;
		idx_tmp = 0;
		num_detected = 0;
		num_detected_tmp = 0;
	}

	public void init(int node_num, boolean is_heu, boolean is_zero) {
		valid = new boolean[node_num*2];
		valid0 = new boolean[node_num*2];
		major = new int[node_num * 2];
		major_tmp = new int[node_num*2];
		four_cycle = new int[4];
		start = new int[3][2];
		tria = new int[3]; // this stores the
		pI = new int[3]; // vertices of distance 1
		pII = new int[3][3];// vertices of distance 2
		pIII = new int[3][3];// vertices of distance 3
		po = new int[3];
		vn = new int[3];
		temp = new int[8];
		left = new int[3];
		right = new int[3];
		idx_major = 0;
		idx_tmp = 0;
		this.as_time = 0;
		this.is_heu = is_heu;
		this.is_zero = is_zero;
	}
}
