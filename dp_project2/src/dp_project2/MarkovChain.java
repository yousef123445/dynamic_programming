package dp_project2;

import java.util.Arrays;

public class MarkovChain {
    private double[][] transitionMatrix;
    private double[] initialVector;
    private int numStates;

    public MarkovChain(double[][] transitionMatrix, double[] initialVector) {
        this.transitionMatrix = transitionMatrix;
        this.initialVector = initialVector;
        this.numStates = transitionMatrix.length;
    }

    public String classifyChain() {
        if (Arrays.stream(transitionMatrix).anyMatch(row -> Arrays.stream(row).anyMatch(p -> p == 1.0))) {
            if (Arrays.stream(transitionMatrix).allMatch(row -> row[Arrays.binarySearch(row, 1.0)] == 1.0)) {
                return "Absorbing";
            } else {
                return "Irreducible, Aperiodic";
            }
        }

        if (checkPeriodicity()) {
            return "Aperiodic";
        }

        if (checkErgodicity()) {
            return "Ergodic";
        }

        return "Unknown";
    }

    private boolean checkPeriodicity() {
        int period = 1;
        double[][] currentPower = transitionMatrix.clone();
        while (!Arrays.deepEquals(currentPower, getIdentityMatrix(numStates))) {
            currentPower = multiplyMatrices(currentPower, transitionMatrix);
            period++;
            if (period > numStates) {
                return false;
            }
        }
        return period % 2 == 0;
    }

    private boolean checkErgodicity() {
        double[][] power = getIdentityMatrix(numStates);
        for (int i = 0; i < numStates; i++) {
            power = multiplyMatrices(power, transitionMatrix);
        }
        return Arrays.stream(power).flatMapToDouble(Arrays::stream).allMatch(p -> p > 0);
    }

    public String classifyState(int stateIndex) {
        if (Arrays.stream(transitionMatrix[stateIndex]).allMatch(p -> p == 0.0) && transitionMatrix[stateIndex][stateIndex] == 1.0) {
            return "Absorbing";
        } else if (checkRecurrentState(stateIndex)) {
            return "Recurrent";
        } else {
            return "Transient";
        }
    }

    private boolean checkRecurrentState(int stateIndex) {
        double[][] power = getIdentityMatrix(numStates);
        for (int i = 0; i < numStates; i++) {
            power = multiplyMatrices(power, transitionMatrix);
        }
        return power[stateIndex][stateIndex] > 0;
    }

    public double[] steadyStateProbability() {
        double[][] transposedMatrix = transposeMatrix(transitionMatrix);
        double[] eigenvalues = new double[numStates];
        double[][] eigenvectors = new double[numStates][numStates];
        computeEigenDecomposition(transposedMatrix, eigenvalues, eigenvectors);

        int indexOfOne = -1;
        for (int i = 0; i < numStates; i++) {
            if (Math.abs(eigenvalues[i] - 1.0) < 1e-10) {
                indexOfOne = i;
                break;
            }
        }

        double[] steadyState = new double[numStates];
        double sum = 0.0;
        for (int i = 0; i < numStates; i++) {
            steadyState[i] = eigenvectors[i][indexOfOne];
            sum += steadyState[i];
        }

        for (int i = 0; i < numStates; i++) {
            steadyState[i] /= sum;
        }

        return steadyState;
    }

    public double[] stateProbabilitiesAfterNSteps(int n) {
        double[][] power = getIdentityMatrix(numStates);
        for (int i = 0; i < n; i++) {
            power = multiplyMatrices(power, transitionMatrix);
        }
        return multiplyMatrixByVector(power, initialVector);
    }

    private double[][] getIdentityMatrix(int size) {
        double[][] identity = new double[size][size];
        for (int i = 0; i < size; i++) {
            identity[i][i] = 1.0;
        }
        return identity;
    }

    private double[][] transposeMatrix(double[][] matrix) {
        double[][] transposed = new double[matrix[0].length][matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                transposed[j][i] = matrix[i][j];
            }
        }
        return transposed;
    }

    private void computeEigenDecomposition(double[][] matrix, double[] eigenvalues, double[][] eigenvectors) {
        // Use a library like JScience or Apache Commons Math to compute the eigendecomposition
        // This is a simplified implementation that may not work for all cases
        for (int i = 0; i < matrix.length; i++) {
            eigenvalues[i] = matrix[i][i];
            for (int j = 0; j < matrix.length; j++) {
                eigenvectors[i][j] = matrix[j][i];
            }
        }
    }

    private double[][] multiplyMatrices(double[][] a, double[][] b) {
        int aRows = a.length;
        int aCols = a[0].length;
        int bCols = b[0].length;

        double[][] result = new double[aRows][bCols];

        for (int i = 0; i < aRows; i++) {
            for (int j = 0; j < bCols; j++) {
                for (int k = 0; k < aCols; k++) {
                    result[i][j] += a[i][k] * b[k][j];
                }
            }
        }

        return result;
    }

    private double[] multiplyMatrixByVector(double[][] matrix, double[] vector) {
        int rows = matrix.length;
        int cols = matrix[0].length;

        double[] result = new double[rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                result[i] += matrix[i][j] * vector[j];
            }
        }

        return result;
    }
    public static void main(String[] args) {
    	
    }
}
