#include <string>
#include <unistd.h>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <math.h>
#include "Eigen/Dense"

#define V_REP_IP_ADDRESS "127.0.0.1"
#define V_REP_PORT 19997//1999;
#define PI 3.14
extern "C" {
#include "extApi.h"
    /*  #include "extApiCustom.h" if you wanna use custom remote API functions! */
}

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

// matrizes para o calculo da distribuição multivariada
MatrixXf grad_f_rl(3,2);
MatrixXf grad_f_p(3,3);
MatrixXf sigma(3,3);

simxInt ddRobotHandle;
simxInt leftMotorHandle;
simxInt rightMotorHandle;
simxInt sensorHandle;
simxInt graphOdometryHandle;
simxFloat pos[3], orientation[3];
int g = 1;

simxFloat goals[8][3] = {{-1.0, 1.5, 3.14/2}, {-0.25, 0.55, 0}, {1, -0.3, -3.14/2}, {1.0, -1.5, -3.14/2}, {1.0, 1.5, 3.14/2}, {0.25, -0.5, -3.14}, {-1.0, -1.5, -3.14/2}, {0.35, 0, 3.14/2}};

simxInt goal= 0;
int state = 0;
simxFloat k_rho = 0.1, k_alpha = 2, k_beta = -0.5;
simxFloat r = 0.02, l = 0.1, b = 0.1, tam_reg = 2;

void getPosition(int clientID, simxFloat pos[]) { //[x,y,theta]

    simxInt ret = simxGetObjectPosition(clientID, ddRobotHandle, -1, pos, simx_opmode_oneshot_wait);
    if (ret > 0) {
        printf("Error reading robot position\n");
        return;
    }

    simxFloat orientation[3];
    ret = simxGetObjectOrientation(clientID, ddRobotHandle, -1, orientation, simx_opmode_oneshot_wait);
    if (ret > 0) {
        printf("Error reading robot orientation\n");
        return;
    }

    simxFloat theta = orientation[2];
    pos[2] = theta;
}

simxInt getSimTimeMs(int clientID) { //In Miliseconds
    return simxGetLastCmdTime(clientID);
}

simxFloat to180range(simxFloat angulo){
    angulo = angulo*180.0/PI;
    if(angulo <= -180){
        angulo = angulo + 360.0;
    }if(angulo >180){
         angulo = angulo - 360.0;
    }
    angulo = (PI*angulo)/180.0;
    return angulo;
}

float to_positive_angle(float angle) {

    angle = fmod(angle, 2 * M_PI);
    while (angle < 0) {
        angle = angle + 2 * M_PI;
    }
    return angle;
}

float smallestAngleDiff(float target, float source) {
    float a;
    a = to_positive_angle(target) - to_positive_angle(source);

    if (a > M_PI) {
        a = a - 2 * M_PI;
    } else if (a < -M_PI) {
        a = a + 2 * M_PI;
    }
    return a;
}

void readOdometers(int clientID, simxFloat &dPhiL, simxFloat &dPhiR) {
    //old joint angle position
    static simxFloat lwprev=0; 
    static simxFloat rwprev=0;
    
    //current joint angle position
    simxFloat lwcur=0;
    simxFloat rwcur=0;

    simxGetJointPosition(clientID, leftMotorHandle, &lwcur, simx_opmode_oneshot);
    simxGetJointPosition(clientID, rightMotorHandle, &rwcur, simx_opmode_oneshot);

    dPhiL = smallestAngleDiff(lwcur, lwprev);
    dPhiR = smallestAngleDiff(rwcur, rwprev);
    lwprev = lwcur;
    rwprev = rwcur;
}

void setTargetSpeed(int clientID, simxFloat phiL, simxFloat phiR) {
    simxSetJointTargetVelocity(clientID, leftMotorHandle, phiL, simx_opmode_oneshot);
    simxSetJointTargetVelocity(clientID, rightMotorHandle, phiR, simx_opmode_oneshot);   
}

inline double to_deg(double radians) {
    return radians * (180.0 / M_PI);
}


void nextGoal(){
    goal = goal+1;
    if (goal>7){
        goal =0;
    }
    state = 0;
}

void motionControl(int clientID, simxFloat &phiL, simxFloat &phiR)
{
    
    simxFloat v, omega, wR, wL;
    simxFloat dx, dy, dtheta, atg, rho;
    simxFloat theta, beta, alpha;
    
    //Localizacao atual do robo [x, y, z]

    getPosition(clientID, pos);
    
    //Angulo atual do robo
    if(simxGetObjectOrientation(clientID, ddRobotHandle, -1, orientation, simx_opmode_oneshot_wait) > 0) {
        printf("Error reading robot orientation\n");
        return;
    }
    
    theta = orientation[2];
    
    //Computar dx, dy, rho, alpha e beta
    dx = goals[goal][0] - pos[0];
    dy = goals[goal][1] - pos[1];
    dtheta = smallestAngleDiff(goals[goal][2], theta);
    rho = sqrt(dx*dx + dy*dy);
    
    atg = atan2(dy, dx);
    atg = to180range(atg);
    alpha = smallestAngleDiff(atg,theta);
    alpha = to180range(alpha);
    beta = (goals[goal][2])-theta - alpha;
    beta = to180range(beta);

    v = k_rho*rho;

    if (v<0.2) v = 0.2;
    
    omega = k_alpha*alpha + k_beta*(beta);

    wR = v + l*omega;
    wL = v - l*omega;

    phiR = wR/r;
    phiL = wL/r;

    // o controle eh instavel com rho ~= 0
    if (rho < 0.05 || state > 0){
        state = 1;
        // faz o robo apenas girar para corrigir theta
        if (abs(dtheta)>0.5){
            phiR = 2*dtheta; // proporcionalmente ao que falta
            phiL = -2*dtheta;
        }
        else{
            phiR = 0;
            phiL = 0;
            state = 2;
        }
    }
}


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

    // cout << endl << "delta_S_R: " << delta_S_R << endl;
    // cout << endl << "delta_S_L: " << delta_S_L << endl;

    MatrixXf cov(2,2);
    cov(0,1) = cov(1,0) = 0;
    cov(0,0) = Kr*abs(delta_S_R);
    cov(1,1) = Kl*abs(delta_S_L);
    return cov;
}

void Cal_sigma(float dPhiL, float dPhiR){
    
    MatrixXf v_cov = get_Cov(dPhiL*r, dPhiR*r);
    sigma = grad_f_p*sigma*(grad_f_p.transpose())+grad_f_rl*v_cov*(grad_f_rl.transpose());

    // cout << endl << "sigma: " << sigma << endl;
    cout << endl << "v_cov: " << v_cov << endl;
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
    MatrixXf x(1,3); 
    x(0,0) = n_x[0]; x(0,1) = n_x[1]; x(0,2) = n_x[2];
    MatrixXf mi(1,3); 
    mi(0,0) = n_mi[0]; mi(0,1) = n_mi[1]; mi(0,2) = n_mi[2];

    //Vector3f x(n_x[0],n_x[1],n_x[2]);
    //Vector3f mi(n_mi[0],n_mi[1],n_mi[2]);

    float y = pow((2*PI),(1.5))*pow((sigma.determinant()),0.5);
    MatrixXf z = (-0.5)*(x-mi)*(sigma.inverse())*((x-mi).transpose());

    cout << "x: " << x << endl;
    cout << "mi: " << mi << endl;
    cout << "sigma: " << sigma << endl;
    cout << "sigma_inv: " << sigma.inverse() << endl;
    cout << "y: " << y << endl;
    cout << "z: " << z(0,0) << endl;
    cout << "exp: " << exp(z(0,0)) << endl;
    return (1.0/y)*exp(z(0,0));
}



void print_reg_belief(int *p){
    for(int i=-tam_reg; i<=tam_reg; i++){
        for(int j=-tam_reg; j<=tam_reg; j++){
            for(int k=-tam_reg; k<=tam_reg; k++){
                cout << belief[p[0]+i][p[1]+j][p[2]+k];
            }
            cout << endl;
        }
        cout << endl;
    }
    cout << endl << endl;
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
    cout << "new belief" << endl;
    cout << pos_belief[0] << " " << pos_belief[1] << " " << pos_belief[2] << endl;
    // nova coordenadas
    float new_coord[3] = {newx,newy,newteta};
    cout << "new_coord: " << new_coord[0] << " " << new_coord[1] << " " << new_coord[2] << endl;

    for(int ni=-tam_reg; ni<=tam_reg; ni++){
        for(int nj=-tam_reg; nj<=tam_reg; nj++){
            for(int nk=-tam_reg; nk<=tam_reg; nk++){
                
                float pos_reg[3] = {newx+ni*TAM_CELL,newy+nj*TAM_CELL,newteta+nk*ANG_MIN};
                cout << "pos_reg: " << pos_reg[0] << " " << pos_reg[1] << " " << pos_reg[2] << endl;
                float bel_reg = 0.0;
                
                for(int i=-tam_reg; i<=tam_reg; i++){
                    for(int j=-tam_reg; j<=tam_reg; j++){
                        for(int k=-tam_reg; k<=tam_reg; k++){
                            bel_reg += dist_mult(pos_reg,new_coord)*belief[x+i][y+j][teta+k];
                            cout << "dist: " << dist_mult(pos_reg,new_coord) << endl;
                            cout << "bel_reg: " << bel_reg << endl;
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

simxFloat Cal_delta_teta(simxFloat dPhiR, simxFloat dPhiL){
    return (dPhiL*r - dPhiR*r)/b;
}


simxFloat Cal_delta_s(simxFloat dPhiR, simxFloat dPhiL){
    return (dPhiL*r +dPhiR*r)/2;
}

simxFloat Cal_deta_x(simxFloat delta_s, simxFloat delta_teta, simxFloat teta){
    return delta_s*cos(teta+(delta_teta/2));
}

simxFloat Cal_deta_y(simxFloat delta_s, simxFloat delta_teta, simxFloat teta){
    return delta_s*sin(teta+(delta_teta/2));
}

int main(int argc, char* argv[]) {

    std::string ipAddr = V_REP_IP_ADDRESS;
    int portNb = V_REP_PORT;

    if (argc > 1) {
        ipAddr = argv[1];
    }

    printf("Iniciando conexao com: %s...\n", ipAddr.c_str());

    int clientID = simxStart((simxChar*) (simxChar*) ipAddr.c_str(), portNb, true, true, 2000, 5);
    
    if (clientID != -1) {
        printf("Conexao efetuada\n");
        
        simxFloat phiL = 0; //rad/s
        simxFloat phiR = 0; //rad/s
        
        //Get handles for robot parts, actuators and sensores:
        simxGetObjectHandle(clientID, "RobotFrame#", &ddRobotHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "LeftMotor#", &leftMotorHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "RightMotor#", &rightMotorHandle, simx_opmode_oneshot_wait);
        simxGetObjectHandle(clientID, "GraphOdometry#", &graphOdometryHandle, simx_opmode_oneshot_wait);

        printf("RobotFrame: %d\n", ddRobotHandle);
        printf("LeftMotor: %d\n", leftMotorHandle);
        printf("RightMotor: %d\n", rightMotorHandle);
        printf("GraphOdometry: %d\n", graphOdometryHandle);

        //start simulation
        int ret = simxStartSimulation(clientID, simx_opmode_oneshot_wait);
        
        if (ret==-1) {
            printf("Não foi possível iniciar a simulação.\n");
            return -1;
        }
        
        printf("Simulação iniciada.\n");

        // posicao inicial
        pos[0] = 1.52; pos[1] = 0.25; pos[2] = 1.09;
        convert_coord_to_ind(pos_belief,pos[0],pos[1],pos[2]);
        belief[pos_belief[0]][pos_belief[1]][pos_belief[2]] = 1.0;

        int count=0;
        //While is connected:
        while (simxGetConnectionId(clientID) != -1) {
            if(state == 2){
                nextGoal();
            }

            //motionControl(clientID, phiL, phiR);
            //setTargetSpeed(clientID, phiL, phiR);
            
            //Read current position:
            //simxFloat pos[3]; //[x,y,theta] in [cm cm rad]
            getPosition(clientID, pos);

            //Read simulation time of the last command:
            simxInt time = getSimTimeMs(clientID); //Simulation time in ms or 0 if sim is not running
            //stop the loop if simulation is has been stopped:
            if (time == 0) break;             
            printf("Posicao: [%.2f %.2f %.2fº], time: %dms\n", pos[0], pos[1],                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      (pos[2]), time);
            
            //Read current wheels angle variation:
            simxFloat dPhiL, dPhiR; //rad
            readOdometers(clientID, dPhiL, dPhiR);
            printf("dPhiL: %.2f dPhiR: %.2f\n", dPhiL*r, dPhiR*r);

            // valores da odometria
            simxFloat delta_teta = Cal_delta_teta(dPhiL,dPhiR);
            simxFloat delta_s = Cal_delta_s(dPhiL,dPhiR);
            simxFloat teta = pos[2]+delta_teta;
            simxFloat delta_x = Cal_deta_x(delta_s,delta_teta,teta);
            simxFloat delta_y = Cal_deta_y(delta_s,delta_teta,teta);

            // erro de odometria
            Cal_grad_f_p(delta_s, teta, delta_teta);
            Cal_grad_f_rl(delta_s, teta, delta_teta);
            if(dPhiL==0.0 && dPhiR==0.0) continue;
            Cal_sigma( dPhiL,  dPhiR);
        
            cout << "pos: " << pos_belief[0] << " " << pos_belief[1] << " " << pos_belief[2] << endl;
            if(count!=0)    print_reg_belief(pos_belief);
            // passo de acao
            action(delta_x, delta_y, delta_teta);
            if(count!=0)    print_reg_belief(pos_belief);

            //Set new target speeds: robot going in a circle:
            //simxFloat phiL = 5; //rad/s
            //simxFloat phiR = 20; //rad/s
            motionControl(clientID, phiL, phiR);
            setTargetSpeed(clientID, phiL, phiR);

            //Let some time for V-REP do its work:
            extApi_sleepMs(2);

            if(count++ > 1){
                //Stop the robot and disconnect from V-Rep;
                setTargetSpeed(clientID, 0, 0);
                simxPauseSimulation(clientID, simx_opmode_oneshot_wait);
                simxFinish(clientID);
                break;
            }
        
        }
        
        //Stop the robot and disconnect from V-Rep;
        setTargetSpeed(clientID, 0, 0);
        simxPauseSimulation(clientID, simx_opmode_oneshot_wait);
        simxFinish(clientID);
        
    } else {
        printf("Nao foi possivel conectar.\n");
        return -2;
    }
    
    return 0;
}


