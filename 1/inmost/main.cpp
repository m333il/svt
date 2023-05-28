#include "inmost.h"
#include <iostream>
#include <chrono>

using namespace INMOST;
        
int main(int argc, char *argv[]) {
	Solver::Initialize(&argc, &argv);

	Sparse::Matrix A;
	std::string matrix_file(argv[1]);
	A.Load(matrix_file);
	std::cout << "Размер системы: " << A.Size() << std::endl;

	Sparse::Vector rhs;
	std::string rhs_file(argv[2]);
	rhs.Load(rhs_file);

	std::string drop_tolerance(argv[3]);

	Solver S(Solver::INNER_ILU2);
	S.SetParameter("verbosity", "2");
	S.SetParameter("drop_tolerance", drop_tolerance);
	
	auto begin = std::chrono::steady_clock::now(); 
	S.SetMatrix(A);
	auto end = std::chrono::steady_clock::now();
	auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "Time for ILU(2): " << elapsed_ms.count() << " ms\n";

	Sparse::Vector sol(rhs);
	begin = std::chrono::steady_clock::now();
	bool solved = S.Solve(rhs, sol);
	end = std::chrono::steady_clock::now();
	elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - begin);
	std::cout << "Time for BiCGStab: " << elapsed_ms.count() << " ms\n";
	std::cout << "Number of iterations: " << S.Iterations() << std::endl;
	std::cout << "Residual: " << S.Residual() << std::endl;
	Solver::Finalize();
	return 0;
}
