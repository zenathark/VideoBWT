package mx.udlap.signalp.wavelet;
import org.opencv.core.*;

public class Wavelet {
	public enum Filters  { db97, db95 };
	
	private Mat coeff;
	private int level;
	private Filters filter;
	
	
	public static void forward(Mat signal, float scale, float[] coeff) {
		int rows = signal.rows();
		int cols = signal.cols();
		Mat result = Mat.zeros(signal.size(), CvType.CV_64F);
		for(int c = 0; c < coeff.length; i+=2) {
			for (int i = 0; i < rows; ++i) {
				for (int j = 1; j < cols; j+=2)
					result.put(i, j, coeff[c] * (signal.get(i,j-1)[0] + signal.get(i,j+1)[0]));
				result.put(i, cols-1, coeff[c] * signal.get(i,cols-1)[0]);
				for (int j = 2; j < cols; j+=2)
					result.put(i, j, coeff[c+1] * (signal.get(i,j-1)[0] + signal.get(i,j+1)[0]));
				result.put(i, 0, coeff[c+1] * signal.get(i,1)[0]);
			}
		}
		for (int i = 0; i < rows; ++i) {
			for (int j = 1; j < cols; j+=2)
				result.put(i, j, scale * result.get(i,j)[0]);
			for (int j = 2; j < cols; j+=2)
				result.put(i, j, result.get(i,j)[0] / scale);
		}
		signal = sortRows(signal);
	}
	
	public Mat sortRows(Mat signal) {
		Mat result = Mat.zeros(signal.size(), signal.type());
		int rows = signal.rows();
		int cols = signal.cols();
		double temp;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols/2; ++j) {
				temp = signal.get(i, j)[0];
			}
		}
	}
	
	public static void inverse(Wavelet wavelet, int level) {
		
	}
}
