import java.math.BigInteger;
import java.io.*;
import java.util.*;
import org.json.simple.*;
import org.json.simple.parser.*;

public class solvePolynomial {

    // Convert string in given base to BigInteger
    static BigInteger convertToDecimal(String value, int base) {
        return new BigInteger(value, base);
    }

    // Solve Vandermonde linear system using Gaussian elimination (with BigDecimal for fractions)
    static BigDecimal[] solveVandermonde(BigDecimal[] xs, BigDecimal[] ys) {
        int n = xs.length;
        BigDecimal[][] A = new BigDecimal[n][n + 1];
        BigDecimal ONE = BigDecimal.ONE;
        BigDecimal ZERO = BigDecimal.ZERO;

        // Build augmented matrix [A|y]
        for (int i = 0; i < n; i++) {
            BigDecimal p = ONE;
            for (int j = n - 1; j >= 0; j--) {
                A[i][j] = p;
                p = p.multiply(xs[i]);
            }
            A[i][n] = ys[i];
        }

        // Gaussian elimination
        for (int i = 0; i < n; i++) {
            // Find pivot
            int pivot = i;
            for (int j = i + 1; j < n; j++) {
                if (A[j][i].abs().compareTo(A[pivot][i].abs()) > 0)
                    pivot = j;
            }
            BigDecimal[] temp = A[i];
            A[i] = A[pivot];
            A[pivot] = temp;

            BigDecimal div = A[i][i];
            for (int k = i; k <= n; k++)
                A[i][k] = A[i][k].divide(div, 50, BigDecimal.ROUND_HALF_UP);

            // Eliminate below
            for (int j = i + 1; j < n; j++) {
                BigDecimal factor = A[j][i];
                for (int k = i; k <= n; k++)
                    A[j][k] = A[j][k].subtract(factor.multiply(A[i][k]));
            }
        }

        // Back substitution
        BigDecimal[] res = new BigDecimal[n];
        for (int i = n - 1; i >= 0; i--) {
            BigDecimal sum = A[i][n];
            for (int j = i + 1; j < n; j++)
                sum = sum.subtract(A[i][j].multiply(res[j]));
            res[i] = sum; // already divided in forward step
        }

        return res; // coefficients [a_n, a_(n-1), ..., c]
    }

    public static void main(String[] args) throws Exception {
        // === Step 1: Load JSON file ===
        JSONParser parser = new JSONParser();
        JSONObject json = (JSONObject) parser.parse(new FileReader("testcase1.json"));

        JSONObject keys = (JSONObject) json.get("keys");
        int n = ((Long) keys.get("n")).intValue();
        int k = ((Long) keys.get("k")).intValue();
        int degree = k - 1;

        // === Step 2: Decode y values ===
        List<BigInteger> yList = new ArrayList<>();
        for (int i = 1; i <= n; i++) {
            JSONObject entry = (JSONObject) json.get(String.valueOf(i));
            if (entry == null) continue;
            int base = Integer.parseInt((String) entry.get("base"));
            String value = (String) entry.get("value");
            yList.add(convertToDecimal(value, base));
        }

        // Select first k points
        BigDecimal[] xs = new BigDecimal[k];
        BigDecimal[] ys = new BigDecimal[k];
        for (int i = 0; i < k; i++) {
            xs[i] = new BigDecimal(i + 1); // x = 1..k
            ys[i] = new BigDecimal(yList.get(i));
        }

        // === Step 3: Solve for coefficients ===
        BigDecimal[] coeffs = solveVandermonde(xs, ys);

        // === Step 4: Extract constant term ===
        BigDecimal constantC = coeffs[coeffs.length - 1];

        System.out.println("Constant term c = " + constantC.toPlainString());
    }
}
