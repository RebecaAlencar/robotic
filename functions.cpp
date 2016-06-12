#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Eigen/Dense"

// medidas do mapa em metros
#define TAM_MAPA 4
#define COMP_PAR_LAT 1
#define COMP_PAR_SUP 2
#define TAM_CELL 0.125
#define ANG_MIN 2.0
#define PI 3.14
// constantes de modelagem do erro
#define Kr 0.2
#define Kl 0.2

using namespace std;
using namespace Eigen;

// modificacao, criando vars
const int ncell = TAM_MAPA/TAM_CELL;
const int gcell = 360/ANG_MIN;

// mapa e belief
bool mapa[ncell][ncell];
float belief[ncell][ncell][gcell];

int pos_belief[3];// indice da posica na matriz belief
float pos[3];// coordenadas reais no mapa

// matrizes para o calculo da distribuição multivariada
MatrixXf grad_f_rl(3,2);
MatrixXf grad_f_p(3,3);
MatrixXf sigma(3,3);

// constantes do robõ
float r = 0.02;
float b = 0.1;

// tamanho da região de busca no belief
int tam_reg = 2;


// monta o mapa
void mount_mapa(){
	// paredes laterais
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

// converte as coordenadas reais do mapa para os indices da matriz do mapa
void convert_coord_to_ind(int *v, float x, float y, float ang){
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

	v[0] = (int)(floor(x/TAM_CELL));
	v[1] = (int)(floor(y/TAM_CELL));
	v[2] = (int)(floor(ang/ANG_MIN));
}


MatrixXf get_Cov(float delta_S_L, float delta_S_R){

	MatrixXf cov(2,2);
	cov(0,1) = cov(1,0) = 0;
	cov(0,0) = Kr*abs(delta_S_R);
	cov(1,1) = Kl*abs(delta_S_L);
	return cov;
}

void Cal_sigma(float dPhiL, float dPhiR){
	
	MatrixXf v_cov = get_Cov(dPhiL*r, dPhiR*r);
	sigma = grad_f_p*sigma*(grad_f_p.transpose())+grad_f_rl*v_cov*(grad_f_rl.transpose());
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

void Cal_grad_f_rl(float delta_s, float teta, float delta_teta){

	grad_f_rl(0,0) = ((1/2)*cos(teta+(delta_teta/2)))-((delta_s/(2*b))*sin(teta+(delta_teta/2)));
	grad_f_rl(0,1) = ((1/2)*cos(teta+(delta_teta/2)))+((delta_s/(2*b))*sin(teta+(delta_teta/2)));
	grad_f_rl(1,0) = ((1/2)*sin(teta+(delta_teta/2)))+((delta_s/(2*b))*cos(teta+(delta_teta/2)));
	grad_f_rl(1,1) = ((1/2)*sin(teta+(delta_teta/2)))-((delta_s/(2*b))*cos(teta+(delta_teta/2)));
	grad_f_rl(2,0) = 1/b;
	grad_f_rl(2,1) = -1/b;

}


float dist_mult(float* n_x, float* n_mi){
	Vector3f x(n_x[0],n_x[1],n_x[2]);
	Vector3f mi(n_mi[0],n_mi[1],n_mi[2]);

	float y = pow((2*PI),(3/2))*pow((sigma.determinant()),1/2);
	float z = (-1/2)*((x-mi).transpose())*(sigma.inverse())*(x-mi);
	return (1/y)*exp(z);
}


void action(float deltax, float deltay, float delta_teta){
	// posicao anterior
	int x = pos_belief[0];
	int y = pos_belief[1];
	int teta = pos_belief[2];

	// posicao depois da odometria
	float newx = pos[0]+deltax;
	float newy = pos[1]+deltay;
	float newteta = pos[2]+delta_teta;
	// nova posica
	convert_coord_to_ind(pos_belief,newx,newy,newteta);
	// nova coordenadas
	float new_coord[3] = {newx,newy,newteta};

	for(int ni=-tam_reg; ni<=tam_reg; ni++){
		for(int nj=-tam_reg; nj<=tam_reg; nj++){
			for(int nk=-tam_reg; nk<=tam_reg; nk++){
				
				float pos_reg[3] = {newx+ni*TAM_CELL,newy+nj*TAM_CELL,newteta+nk*ANG_MIN};
				float bel_reg = 0.0;
				
				for(int i=-tam_reg; i<=tam_reg; i++){
					for(int j=-tam_reg; j<=tam_reg; j++){
						for(int k=-tam_reg; k<=tam_reg; k++){
							bel_reg += dist_mult(pos_reg,new_coord)*belief[x+i][y+j][teta+k];
						}
					}
				}
				
				int pos_aux[3];	
				convert_coord_to_ind(pos_aux,pos_reg[0],pos_reg[1],pos_reg[2]);
				belief[pos_aux[0]][pos_aux[1]][pos_aux[2]] = bel_reg;
			}
		}
	}
}


int main(){
	// posição inicial do robo
	float x = 1.5, y = 0.25, ang = 0.0;
	// constroi o mapa
	mount_mapa();
	// belief inicial
	convert_coord_to_ind(pos_belief,x,y,ang);
	belief[pos_belief[0]][pos_belief[1]][pos_belief[2]] = 1.0;

}	






