#include <iomanip>
#include <fstream>
#include <iostream>
#include <vector>
#include <list>

using namespace std;

struct Cell
{
	double c;
	double x;
	bool fromBasis;

	Cell()
	{
		c = 0;
		x = 0;
		fromBasis = false;
	}
};

double sum(vector<double>& vector)
{
	double sum = 0;

	for(double item : vector)
	{
		sum += item;
	}

	return sum;
}

void setBalance(vector<vector<Cell>>& transportMatrix, vector<double>& amounts1, vector<double>& amounts2)
{
	double sum1 = sum(amounts1);
	double sum2 = sum(amounts2);
	
	if(sum1 > sum2)
	{
		for(vector<Cell>& vector : transportMatrix)
		{
			vector.resize(vector.size() + 1);
		}

		amounts2.resize(amounts2.size() + 1, sum1 - sum2);
	}
	else
	if(sum1 < sum2)
	{
		transportMatrix.resize(transportMatrix.size() + 1);
		transportMatrix[transportMatrix.size() - 1].resize(transportMatrix[0].size());
		amounts1.resize(amounts1.size() + 1, sum2 - sum1);
	}
}

void createFirstPlan(vector<vector<Cell>>& transportMatrix, vector<double>& amounts1, vector<double>& amounts2, vector<vector<int>>& basis)
{
	int latitude;
	int longtitude;
	double delivery = 0;
	int counterOfCells = 0;
	int counterOfCellsInBasis = 0;
	vector<vector<Cell>> transportMatrixTemp(transportMatrix.size());
	vector<double> amounts1Temp;
	vector<double> amounts2Temp;

	list<int> rows;
	list<int> columns;

	while(counterOfCellsInBasis < transportMatrix.size() + transportMatrix[0].size() - 1)
	{
		for(int i = 0; i < transportMatrix.size(); i++)
		{
			transportMatrixTemp[i] = transportMatrix[i];
		}

		amounts1Temp = amounts1;
		amounts2Temp = amounts2;
	
		for(int i = 0; i < transportMatrix.size(); i++)
		{
			rows.push_back(i);
		}

		for(int j = 0; j < transportMatrix[0].size(); j++)
		{
			columns.push_back(j);
		}

		for(vector<int>& vector : basis)
		{
			vector.resize(0);
		}
		
		latitude	= counterOfCells / transportMatrix[0].size();
		longtitude  = counterOfCells % transportMatrix[0].size();
		counterOfCellsInBasis = 0;
		
		while(true)
		{
			counterOfCellsInBasis++;
			basis[latitude].push_back(longtitude);
			delivery = min(amounts1Temp[latitude], amounts2Temp[longtitude]);
			amounts1Temp[latitude]	 -= delivery;
			amounts2Temp[longtitude] -= delivery;
			transportMatrixTemp[latitude][longtitude].x			= delivery;
			transportMatrixTemp[latitude][longtitude].fromBasis = true;

			if(amounts1Temp[latitude] == 0)
			{
				rows.erase(find(rows.begin(), rows.end(), latitude));
			}
			
			if(amounts2Temp[longtitude] == 0)
			{
				columns.erase(find(columns.begin(), columns.end(), longtitude));
			}

			if(!rows.empty())
			{
				latitude = rows.front();
				
				if(!columns.empty())
				{
					longtitude = columns.front();
				}
			}
			else
			{
				if(!columns.empty())
				{
					longtitude = columns.front();
				}
				else
				{
					break;
				}
			}
		}

		counterOfCells++;
	}

	for(int i = 0; i < transportMatrixTemp.size(); i++)
	{
			transportMatrix[i] = transportMatrixTemp[i];
	}
}

void calculatePotentialsForLines(vector<vector<Cell>>& transportMatrix, vector<bool>& isCalculated, vector<double>& u, vector<double>& v, vector<vector<int>>& basis, int num)
{
	for(int j : basis[num])
	{
		v[j] = -transportMatrix[num][j].c - u[num];
	}

	isCalculated[num] = true;

	for(int j : basis[num])
	{
		for(int i = 0; i < basis.size(); i++)
		{
			if(find(basis[i].begin(), basis[i].end(), j) != basis[i].end() && isCalculated[i] == false)
			{
				u[i] = -transportMatrix[i][j].c - v[j];
				calculatePotentialsForLines(transportMatrix, isCalculated, u, v, basis, i);
			}
		}
	}
}

void calculatePotentials(vector<vector<Cell>>& transportMatrix, vector<double>& u, vector<double>& v, vector<vector<int>>& basis)
{
	vector<bool> isCalculated(u.size());

	calculatePotentialsForLines(transportMatrix, isCalculated, u, v, basis, 0);
}

bool checkEstimations(vector<vector<Cell>>& transportMatrix, vector<double>& u, vector<double>& v, pair<int, int>& cellOfMaxEstimation)
{
	double estimation;
	double maxEstimation = 0;
	
	for(int i = 0; i < transportMatrix.size(); i++)
	{
		for(int j = 0; j < transportMatrix[0].size(); j++)
		{
			if(transportMatrix[i][j].fromBasis == false)
			{
				estimation = -transportMatrix[i][j].c - (u[i] + v[j]);

				if(estimation > maxEstimation)
				{
					maxEstimation = estimation;
					cellOfMaxEstimation.first	= i;
					cellOfMaxEstimation.second	= j;
				}
			}
		}
	}

	if(maxEstimation == 0)
	{
		return true;
	}
	else
	{
		return false;
	}
}

bool findCycle(vector<vector<Cell>>& transportMatrix, pair<int, int>& startCell, pair<int, int>& cellOfMaxEstimation, vector<pair<int, int>>& cycle, bool condition)
{
	pair<int, int> currCell;
	int limit;

	currCell.first	= startCell.first;
	currCell.second = startCell.second;
	int* variableParameter;

	if(condition)
	{
		limit = transportMatrix.size();
		variableParameter = &currCell.first;
	}
	else
	{
		limit = transportMatrix[0].size();
		variableParameter = &currCell.second;
	}

	for(int i = 1; i < limit; i++)
	{
		*variableParameter = (++(*variableParameter)) % limit;

		if(currCell == cellOfMaxEstimation)
		{
			return true;
		}
		
		if(transportMatrix[currCell.first][currCell.second].fromBasis == true)
		{
			if(findCycle(transportMatrix, currCell, cellOfMaxEstimation, cycle, !condition))
			{
				cycle.push_back(currCell);
				return true;
			}
		}
	}

	return false;
}

void findMinimumOfAmounts(vector<vector<Cell>>& transportMatrix, vector<pair<int, int>>& cycle, pair<int, int>& cellOfMinAmount)
{
	double minAmount = transportMatrix[cycle[0].first][cycle[0].second].x;
	cellOfMinAmount = cycle[0];

	for(int i = 2; i < cycle.size(); i += 2)
	{
		if(minAmount > transportMatrix[cycle[i].first][cycle[i].second].x)
		{
			minAmount = transportMatrix[cycle[i].first][cycle[i].second].x;
			cellOfMinAmount = cycle[i];
		}
	}
}

void doModifications(vector<vector<Cell>>& transportMatrix, vector<pair<int, int>>& cycle, vector<vector<int>>& basis, pair<int, int>& cellOfMaxEstimation, pair<int, int>& cellOfMinAmount)
{
	basis[cellOfMinAmount.first].erase(find(basis[cellOfMinAmount.first].begin(), basis[cellOfMinAmount.first].end(), cellOfMinAmount.second));
	basis[cellOfMaxEstimation.first].push_back(cellOfMaxEstimation.second);
	transportMatrix[cellOfMinAmount.first][cellOfMinAmount.second].fromBasis = false;
	transportMatrix[cellOfMaxEstimation.first][cellOfMaxEstimation.second].fromBasis = true;

	double increment = transportMatrix[cellOfMinAmount.first][cellOfMinAmount.second].x;

	transportMatrix[cellOfMaxEstimation.first][cellOfMaxEstimation.second].x += increment;

	for(int i = 0; i < cycle.size(); i++)
	{
		transportMatrix[cycle[i].first][cycle[i].second].x += pow(-1, i + 1) * increment;
	}
}

void clean(vector<pair<int, int>>& cycle, vector<double>& u)
{
	u[0] = 0;
	cycle.resize(0);
}

void fitMatrix(vector<vector<Cell>>& transportMatrix, int n1, int n2)
{
	transportMatrix.resize(n1);
	for(vector<Cell>& vector : transportMatrix)
	{
		vector.resize(n2);
	}
}

int main()
{
	int n1;
	int n2;
	
	ifstream fin("input.txt");
	fin >> n1 >> n2;

	vector<vector<Cell>> transportMatrix(n1);
	vector<double> amounts1(n1);
	vector<double> amounts2(n2);
	vector<vector<int>> basis;
	vector<double> u;
	vector<double> v;
	pair<int, int> cellOfMaxEstimation;
	vector<pair<int, int>> cycle;
	pair<int, int> cellOfMinAmount;

	for(vector<Cell>& vector : transportMatrix)
	{
		vector.resize(n2);
		for(Cell& item : vector)
		{
			fin >> item.c;
		}
	}

	for(double& item : amounts1)
	{
		fin >> item;
	}

	for(double& item : amounts2)
	{
		fin >> item;
	}

	fin.close();

	setBalance(transportMatrix, amounts1, amounts2);
	
	basis.resize(amounts1.size());
	u.resize(amounts1.size());
	v.resize(amounts2.size());
	
	createFirstPlan(transportMatrix, amounts1, amounts2, basis);
	calculatePotentials(transportMatrix, u, v, basis);

	while(!checkEstimations(transportMatrix, u, v, cellOfMaxEstimation))
	{
		findCycle(transportMatrix, cellOfMaxEstimation, cellOfMaxEstimation, cycle, true);
		findMinimumOfAmounts(transportMatrix, cycle, cellOfMinAmount);
		doModifications(transportMatrix, cycle, basis, cellOfMaxEstimation, cellOfMinAmount);
		clean(cycle, u);
		calculatePotentials(transportMatrix, u, v, basis);
	}

	fitMatrix(transportMatrix, n1, n2);

	ofstream fout("output.txt");

	for(vector<Cell>& vector : transportMatrix)
	{
		for(Cell& item : vector)
		{
			fout << fixed << setprecision(3) << setw(20) << item.x;
		}
		fout << "\n";
	}

	return 0;
}