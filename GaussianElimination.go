package main

import (
	"errors"
	"fmt"
	"math"
)

type Matrix [][]float64

func main() {

	var a = [][]float64{
		{1.0, 2.0, 3.0, 1.0},
		{4.0, 5.0, 6.0, 1.0},
		{1.0, 0.0, 1.0, 1.0},
	}
	fmt.Println(ComputeP(a))

}

//CalculateP takes in D and L matrix, P2=0. Solve the remaining n-1 P by n equations. Returning p matrix is a list of float64 variables.
func ComputeP(A Matrix) ([]float64, error) {
	N := len(A)
	//solve linear equations by guassian GaussianElimination : AX=b
	//Initialize the faction matrix dij/lij, augment A by adding the b terms.

	for i := range A {
		//Find pivot for column i
		iMax := i
		max := math.Abs(A[i][i])
		for j := i + 1; j <= N-1; j++ {
			if temp := math.Abs(A[j][i]); temp > max {
				iMax = j
				max = temp
			}
		}
		if A[iMax][i] == 0 {
			return nil, errors.New("singular")
		}
		//swap row A[iMax] and A[i] Or A[i],A[imax]=A[imax],A[i]
		tempr := A[i]
		A[i] = A[iMax]
		A[iMax] = tempr

		//for all rows below pivot
		for r := i + 1; r <= N-1; r++ {
			//Fill 0 with the lower part of the pivot column
			A[r][i] = 0

			f := -A[r][i] / A[i][i]
			//for all remaining nonzero elements in current row
			for c := i + 1; c <= N; c++ {
				A[r][c] += f * A[i][c]
			}
		}
	}

	//retrive p from the upper triangular matrix
	p := make([]float64, N)
	for i := N - 1; i >= 0; i-- {
		p[i] = A[i][N]
		for j := i + 1; j <= N-1; j++ {
			p[i] -= A[i][j] * p[j]
		}
		p[i] /= A[i][i]
	}

	return p, nil
}
