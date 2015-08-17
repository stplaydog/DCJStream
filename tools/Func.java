package tools;

import java.text.DecimalFormat;

public class Func {

	private static final int ASCII_VALUE_OF_ZERO = 48;

	public static int min(int start[][]) {
		if (start[0][1] <= start[1][1]) {
			if (start[0][1] <= start[2][1])
				return 0;
			else
				return 2;
		} else if (start[1][1] <= start[2][1])
			return 1;
		else
			return 2;
	}

	public static char[] itoa(int number) {
		boolean negative = false;
		if (number < 0) {
			negative = true;
			number = 0 - number;
		}

		if (number >= 0 && number <= 9) {
			char temp = (char) (ASCII_VALUE_OF_ZERO + number);
			if (!negative) {
				return new char[] { temp };
			}

			return new char[] { '-', temp };
		}

		// define an array of which can hold 12 characters
		// the max integer is 10 digits long - 1 for negative character
		char[] digits = new char[12];

		// now let's divide the number by 10 and keep adding the remainder
		int digitPosition = 0;

		do {
			int remainder = number % 10;
			number = number / 10;

			digits[digitPosition++] = (char) (ASCII_VALUE_OF_ZERO + remainder);
		} while (number > 0);

		// add negative sign if needed
		if (negative) {
			digits[digitPosition++] = '-';
		}

		// now reverse the array
		for (int i = 0; i < digitPosition / 2; i++) {
			char temp = digits[i];
			digits[i] = digits[digitPosition - i - 1];
			digits[digitPosition - i - 1] = temp;
		}

		return digits;
	}

	public static double to_sec(long timemillis) {
		return (double) timemillis / (double) 1000;
	}

	public static double rt(double d) {
		DecimalFormat twoDForm = new DecimalFormat("#.##");
		return Double.valueOf(twoDForm.format(d));
	}
}
