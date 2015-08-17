package order;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;

public class GeneOrderCir extends GeneOrder {

	public void init(String file) {
		boolean init = false;
		// TODO Auto-generated constructor stub
		name = new String[3];
		int gene_id = 0;
		try {
			BufferedReader reader = new BufferedReader(new FileReader(new File(
					file)));
			String line = reader.readLine();
			while (line != null) {
				line = line.trim();
				if (line.indexOf(">") != -1)
					name[gene_id] = line.replace(">", "");
				else {
					String genes[] = line.split(" ");
					if (init == false) {
						init = true;
						gene_order = new int[genes.length][3];
						this.gene_num = genes.length;
					}
					for (int i = 0; i < genes.length; i++) {
						gene_order[i][gene_id] = Integer.parseInt(genes[i]);
					}
					gene_id++;
				}
				line = reader.readLine();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public void init(int gene_num, int chr_num) {
		this.gene_num = gene_num;
		gene_order = new int[gene_num][3];
		name = new String[3];

		for (int i = 0; i < 3; i++) {
			name[i] = String.valueOf(i);
			for (int j = 0; j < gene_num; j++) {
				gene_order[j][i] = (j + 1);
			}
		}
	}

	public String toString() {
		StringBuffer result = new StringBuffer();
		for (int i = 0; i < 3; i++) {
			result.append(">" + this.name[i] + "\n");
			for (int j = 0; j < this.gene_num; j++)
				result.append(this.gene_order[j][i] + " ");
			result.append("\n");
		}
		return result.toString();
	}

	@Override
	public String toString(int genome_id) {
		// TODO Auto-generated method stub
		StringBuffer result = new StringBuffer();
		for (int i = 0; i < this.gene_num; i++)
			result.append(this.gene_order[i][genome_id] + " ");
		return result.toString().trim();
	}

}
