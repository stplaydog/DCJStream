package order;

public abstract class GeneOrder {
	public int gene_order[][];
	public int gene_num;
	public int chr_num[];
	public String name[];

	public abstract String toString();

	public abstract void init(String file);
	
	public abstract void init(int gene_num, int chr_num);
	
	public abstract String toString(int genome_id);
}
