#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <cmath>

#define tam_mapa 4
#define comp_parede_lat 1
#define comp_parede_sup 2
#define tam_cell 0.125
#define disc_angle 2.0

#define Kr 0.2
#define Kl 0.2

using namespace std;

// mapa
bool mapa[32][32];
double belief[32][32][int(360/disc_angle)];
double grad_f_rl[3][2];
double grad_f_p[3][3];
double sigma[3][3];
// monta o mapa
void mount_mapa(){
	// paredes laterais
	int ncell = 32;
	for(int i=0; i<ncell; i++){
		mapa[0][i] = true;
		mapa[ncell-1][i] = true;
		mapa[i][0] = true;
		mapa[i][ncell-1] = true;
	}

	// paredes internas menores
	int comp_par = comp_parede_lat/tam_cell;
	int pos_mid = (tam_mapa/tam_cell)/2;
	for(int i=0; i<=comp_par; i++){
		mapa[pos_mid][i] = true;
		mapa[pos_mid][ncell-1-i] = true;
	}

	// parede interna maior
	comp_par = comp_parede_sup/tam_cell;
	for(int i=0; i<=comp_par; i++){
		mapa[i][pos_mid] = true;
	}
}

void print_mapa(){
	int ncell = 32;
	for(int i=0; i<ncell; i++){
		cout << i << " ";
		for(int j=0; j<ncell; j++)
			cout << mapa[i][j] << " ";
		cout << endl;
	}
}

// converte para os indices da matriz do mapa
void convert_coord_to_mapa(float *x, float *y){
	if(*x < 0.0 && *y < 0.0 || *x>0.0 && *y>0.0){
		*x = *x+2.0;
		*y = *y+2.0;
	}else if(*x<0.0 && *y>0.0){
		*x = 2.0+*x;
		*y = *y+2.0;
	}else if(*x>0.0 && *y<0.0){
		*y = -*y+2.0;
	}
}

void convert_angle(float *ang){
	*ang = floor(*ang/disc_angle);
}


void Cal_grad_f_p(float delta_s, float teta, float delta_teta){

	grad_f_p[0][0]= 1;
	grad_f_p[0][1]= 0;
	grad_f_p[0][2]= -delta_s*sin(teta+(delta_teta/2));
	grad_f_p[1][0]= 0;
	grad_f_p[1][1]= 1;
	grad_f_p[1][2]= delta_s*sin(teta+(delta_teta/2));
	grad_f_p[2][0]= 0;
	grad_f_p[2][1]= 0;
	grad_f_p[2][2]= 1;

}



void Cal_sigma(){

	sigma = 

}

void Cal_grad_f_rl(float delta_s, float teta, float delta_teta, float b){

	grad_f_rl[0][0] = ((1/2)*cos(teta+(delta_teta/2)))-((delta_s/(2*b))*sin(teta+(delta_teta/2)));
	grad_f_rl[0][1] = ((1/2)*cos(teta+(delta_teta/2)))+((delta_s/(2*b))*sin(teta+(delta_teta/2)));
	grad_f_rl[1][0] = ((1/2)*sin(teta+(delta_teta/2)))+((delta_s/(2*b))*cos(teta+(delta_teta/2)));
	grad_f_rl[1][1] = ((1/2)*sin(teta+(delta_teta/2)))-((delta_s/(2*b))*cos(teta+(delta_teta/2)));
	grad_f_rl[2][0] = 1/b;
	grad_f_rl[2][1] = -1/b;



}

int main(){
	// constroi o mapa
	mount_mapa();

	// posição do robo
	float x = 1.5, y = 0.25, ang = 0.0;
	convert_coord_to_mapa(&x,&y);
	convert_angle(&ang);
	// inicial
	int indx = (int)(x/tam_cell);
	int indy = (int)(y/tam_cell);
	int indang = (int)(ang);
	mapa[indx][indy] = true;
	belief[indx][indy][indang] = 1.0;

	print_mapa();
}