#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include <Eigen/Dense>

// medidas do mapa em metros
#define TAM_MAPA 4
#define COMP_PAR_LAT 1
#define COMP_PAR_SUP 2
#define TAM_CELL 0.125
#define ANG_MIN 2.0

// constantes de modelagem do erro
#define Kr 0.2
#define Kl 0.2

using namespace std;
using namespace Eigen;

// modificacao, criando vars
const int ncell = TAM_MAPA/TAM_CELL;
const int gcell = 360/ANG_MIN;

// mapa
bool mapa[ncell][ncell];
double belief[ncell][ncell][gcell];
VectorXf coord_robo_belief(3);
VectorXd coord_robo_mapa(3);
MatrixXf grad_f_rl(3,2);
MatrixXf grad_f_p(3,3);
MatrixXf sigma(3,3);

// monta o mapa
void mount_mapa(){
	// paredes laterais
	int ncell = (int)(TAM_MAPA/TAM_CELL);
	for(int i=0; i<ncell; i++){
		mapa[0][i] = true;
		mapa[ncell-1][i] = true;
		mapa[i][0] = true;
		mapa[i][ncell-1] = true;
	}

	// paredes internas menores
	int comp_par = COMP_PAR_LAT/TAM_CELL;
	int pos_mid = (TAM_MAPA/TAM_CELL)/2;
	for(int i=0; i<=comp_par; i++){
		mapa[pos_mid][i] = true;
		mapa[pos_mid][ncell-1-i] = true;
	}

	// parede interna maior
	comp_par = COMP_PAR_SUP/TAM_CELL;
	for(int i=0; i<=comp_par; i++){
		mapa[i][pos_mid] = true;
	}
}

void print_mapa(){
	for(int i=0; i<ncell; i++){
		cout << i << " ";
		for(int j=0; j<ncell; j++)
			cout << mapa[i][j] << " ";
		cout << endl;
	}
}

// converte para os indices da matriz do mapa
VectorXd convert_coord_to_ind(float x, float y, float ang){
	// x,y
	if(x < 0.0 && y < 0.0 || x>0.0 && y>0.0){
		x = x+2.0;
		y = y+2.0;
	}else if(x<0.0 && y>0.0){
		x = 2.0+x;
		y = y+2.0;
	}else if(x>0.0 && y<0.0){
		y = -y+2.0;
	}

	if(ang < 0.0){
		ang = 360.0 + ang;
	}

	VectorXd v(3);
	v(0) = (int)(floor(x/TAM_CELL));
	v(1) = (int)(floor(y/TAM_CELL));
	v(2) = (int)(floor(ang/ANG_MIN));
	return v;
}

void Cal_grad_f_p(float delta_s, float teta, float delta_teta){

	grad_f_p(0,0)= 1;
	grad_f_p(0,1)= 0;
	grad_f_p(0,2)= -delta_s*sin(teta+(delta_teta/2));
	grad_f_p(1,0)= 0;
	grad_f_p(1,1)= 1;
	grad_f_p(1,2)= delta_s*sin(teta+(delta_teta/2));
	grad_f_p(2,0)= 0;
	grad_f_p(2,1)= 0;
	grad_f_p(2,2)= 1;

}



void Cal_sigma(){
	
}

void Cal_grad_f_rl(float delta_s, float teta, float delta_teta, float b){

	grad_f_rl(0,0) = ((1/2)*cos(teta+(delta_teta/2)))-((delta_s/(2*b))*sin(teta+(delta_teta/2)));
	grad_f_rl(0,1) = ((1/2)*cos(teta+(delta_teta/2)))+((delta_s/(2*b))*sin(teta+(delta_teta/2)));
	grad_f_rl(1,0) = ((1/2)*sin(teta+(delta_teta/2)))+((delta_s/(2*b))*cos(teta+(delta_teta/2)));
	grad_f_rl(1,1) = ((1/2)*sin(teta+(delta_teta/2)))-((delta_s/(2*b))*cos(teta+(delta_teta/2)));
	grad_f_rl(2,0) = 1/b;
	grad_f_rl(2,1) = -1/b;

}

int main(){
	// constroi o mapa
	mount_mapa();
	// posição inicial do robo
	float x = 1.5, y = 0.25, ang = 0.0;
	// belief inicial
	coord_robo_mapa = convert_coord_to_ind(x,y,ang);
	belief[(int)coord_robo_mapa(0)][(int)coord_robo_mapa(1)][(int)coord_robo_mapa(2)] = 1.0;
}






