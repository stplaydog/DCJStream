import java.io.File;

import tools.ResultFilter;

public class Stats {
	public static void main(String args[]) {
		String folder = args[0];
		String type = args[1];

		File fold = new File(folder);
		String files[] = fold.list();

		for (int i = 0; i < files.length; i++) {
			if (type.equals("complexity")) {
				ResultFilter.filter_complexity(files[i], folder);
			} else if (type.equals("complete"))
				ResultFilter.filter_complete(files[i], folder);
			else if (type.equals("parallel"))
				ResultFilter.filter_complete_p(files[i], folder);
		}
	}
}
