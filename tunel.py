import matplotlib.pyplot as plt
import numpy as np

#     x1       _x2________________________________
#             /        o
#            /         o
#           /theta     o
# -----------          o
#                      o
#                      o
# ----------           o
#          |\          o
#          | \         o
#          |  \___________________________________
#-1-      -2-       -3-  -4-                    -5-

# GRAVIDADE

g = 9.81;
pre_1 = 101325;

# COOLER INFO

etta = 1; #Rendimento do cooler
Vaz_cfm = 40.9; #Potencia do Cooler em Cfm

# FLUIDO INFO

rho = 1.169; # densidade do ar, 20 graus C e 1 atm 
mu = 18.4 * 10**(-6) # viscosidade dinamica do ar

# TUNEL INFO

x1 = 0.095; #Diametro hidraulico de entrada do tunel
x2 = 0.5;  #Diametro hidraulico da secao de teste

angle = 45; # Angulo de abertura do difusor

def deg_to_rad(x): # Funcao que converte graus para radianos
    return x * np.pi/180;

###############
    
    
Vaz = 0.0004719474 * Vaz_cfm * etta; #Conversao da vazao de cfm para m3/s

A1 = (x1**2); #Area da entrada do tunel
A2 = (x2**2); #Area da secao de testes

AR = A2/A1; #Razao de area

# VELOCIDADE EM 1

v1 = Vaz / A1;

display(v1);

Re_1 = rho * x1 * v1 / mu;

def project(x1, x2, etta, angle): ###### Esta funcao corresponde ao projeto do tunel
    
    ##############
    
    ###### Funcoes recorrentes durante o projeto
    
    ##############
    
    
    def L_f(v_in, dist, diameter): #funcao para calcular a velocidade apeos a perda de carga por friccao em uma secao reta
        dev = 1;
        f_old = 1;
        
        Re_f = rho * diameter * v_in / mu;
        
        if Re_f < 2100:
            f_l = 64/Re_f;
        else:
            while dev > (10**(-6)):
                f_l = (2 * np.log10(Re_f * np.sqrt(f_old)) - 0.8)**(-2);
                dev = np.abs((f_l - f_old)/f_old);
                f_old = f_l;
             
        h_f = f_l * dist * v_in**2/(diameter * 2); #perda de carga por friccao
        return np.sqrt(2*(v_in**2/2 - h_f)); #velocidade media no ponto 'dist' do tunel\
    
    #########################
    
    def deltaP_f(v_in, dist, diameter): #funcao para calcular diferenca de pressao em um duto reto ocasionado pela perda de carga
        dev = 1;
        f_old = 1;
        
        Re_f = rho * diameter * v_in / mu;
        
        if Re_f < 2100:
            f_l = 64/Re_f;
        else:
            while dev > (10**(-6)):
                f_l = (2 * np.log10(Re_f * np.sqrt(f_old)) - 0.8)**(-2);
                dev = np.abs((f_l - f_old)/f_old);
                f_old = f_l;
             
        return f_l * dist * (v_in**2) * rho / (diameter * 2); #diferenca de pressao no ponto 'dist' do tunel

    ##############
    
    ###### DIMENSIONAMENTO DO CORREDOR INICIAL
    
    ##############
    
    init_dist = 0.10; #comprimento do corredor inicial
    
    dist_vec_init = np.arange(0.01, init_dist, 0.01); #discretizacao do comprimento do corredor inicial
    
    v_vec_init = np.zeros(np.size(dist_vec_init)); #criacao do vetor de velocidade ao longo da secao de teste
    Re_vec_init = np.zeros(np.size(dist_vec_init)); #criacao do vetor de numero de reynolds ao longo da secao de teste
    
    for i in range (0, np.size(dist_vec_init), 1):
        v_vec_init[i] = L_f(v1, dist_vec_init[i], x1); #calculo da velocidade em cada ponto discretizado da secao de teste
        Re_vec_init[i] = rho * x1 * v_vec_init[i] / mu; #calculo do numero de reynolds em cada ponto discretizado da secao de teste
    
    ##############
    
    ###### DIMENSIONAMENTO DO DIFUSOR
    
    ##############
    
    # FATOR DE FRICCAO

    v2 = v_vec_init[np.size(v_vec_init) - 1];
    
    pre_2 = pre_1 - deltaP_f(v1, init_dist, x1); #pressao na entrada do difusor
    
    display(pre_2);
    
    Re_2 = rho * x1 * v2 / mu;

    def friction(): #Funcao para calcuar a perda de carga por friccao no difusor
        dev = 1; #estimativa inicial do desvio relativo para o criterio de parada
        f_old = 1; # estimativa inicial do fator de friccao
        
        if Re_2 < 2100:
            f = 64/Re_2;
        else:
            while dev > (10**(-6)):
                f = (2 * np.log10(Re_2 * np.sqrt(f_old)) - 0.8)**(-2);
                dev = np.abs((f - f_old)/f_old);
                f_old = f;
        
        return f;
    
    # PERDA DE CARGA NO DIFUSOR
    
    f = friction(); # fator de friccao no difusor
    
    Kf = (1 - 1/(AR**2)) * f/(8 * np.sin(deg_to_rad(angle))); # coeficiente de perda de carga por friccao no difusor
    
    def function_Ke(x): #Funcao para calculo do coeficiente de perda de carga angular no difusor
        if x <= 1.5:
            Ke = 0.09623 - 0.004152*x;
        elif x > 1.5 and x < 5:
            Ke = 0.1222 - 0.04590*x + 0.02203*x**2 + 0.0032690*x**3 - 0.0006145*x**4 - 0.00002800*x**5 + 0.00002337*x**6;
        elif x >= 5:
            Ke = -0.01322 + 0.05866 * x;
        return Ke;
    
    Ke = function_Ke(angle); 
    Kex = Ke*(1 - 1/AR)**2; #coefiente de perda de carga angular no difusor
    
    Kd = Kf + Kex; # coeficiente de perda de carga total no difusor
    
    pre_3 = pre_2 - Kd * rho * (v2**2)/2; #pressao na saida do difusor
    
    h_d = Kd * (v2**2)/(2*g); # Perda de carga no difusor
    
    ##############
    
    ###### DIMENSIONAMENTO DO COMPRIMENTO DA COLMEIA
    
    ##############
    
    
    #Perda de carga na colmeia
    
    v3 = A1/A2 * np.sqrt(2*(v2**2/2 - h_d)); #velocidade de entrada na colmeia (saida do difusor)
    
    Re_3 = rho * x2 * v3 / mu; #numero de reynolds na entrada da colmeia
    
    delta = 0.00006; #rugosidade do material da colmeia (canudo de plastico)
    d_h = 0.005; #diametro do canudo
    l_h = 0.10; #comprimento do canudo
    beta_h = 0.8; #coeficiente da porosidade da celula
    
    def function_Kh(v): #funcao para calcular o coeficiente de perda de carga na colmeia
    
        Re_delta = rho * v * d_h / mu; #numero de reynold na entrada de uma celula
    
        if (Re_delta < 275):
            delta_h = 0.375 * ((delta/d_h)**(0.4)) * (Re_delta**(-0.1));
        else:
            delta_h = 0.214 * (delta/d_h)**(0.4);
    
        return delta_h * (l_h/d_h + 3)*(1/beta_h)**2 + (1/beta_h - 1)**2; #coeficiente de perda de carga na colmeia
    
    K_h = function_Kh(v3); #coeficiente de perda de carga na colmeia
    
    h_h = K_h * (v3**2)/(2 * g); #perda de carga na colmeia
    
    pre_4 = pre_3 - K_h * rho * (v2**2)/2;
    
    ##############
    
    ###### DIMENSIONAMENTO DA SECAO DE TESTE
    
    ##############
        
    sec_dist = 0.5; #comprimento da secao de teste
    
    v4 = np.sqrt(2*(v3**2/2 - h_h)); #velocidade na entrada da secao de teste
    Re_4 = rho * x2 * v4 / mu; #numero de Reynolds na entrada da secao de teste
    
    dist_vec = np.arange(0.01, sec_dist, 0.01); #discretizacao do comprimento da secao de teste
    
    v_vec = np.zeros(np.size(dist_vec)); #criacao do vetor de velocidade ao longo da secao de teste
    Re_vec = np.zeros(np.size(dist_vec)); #criacao do vetor de numero de reynolds ao longo da secao de teste
    
    for i in range (0, 49, 1):
        v_vec[i] = L_f(v4, dist_vec[i], x2); #calculo da velocidade em cada ponto discretizado da secao de teste
        Re_vec[i] = rho * x2 * v_vec[i] / mu; #calculo do numero de reynolds em cada ponto discretizado da secao de teste
    
    pre_5 = pre_4 - deltaP_f(v4, sec_dist, x2);
    
    dif_dist = 0.5*(x2-x1)*np.tan(deg_to_rad(angle));
    
    #####PLOT DOS GRAFICOS
        
    fig = plt.figure()
    
    Re_chart = plt.figure(); #Grafico do numero de Reynolds

    plt.plot(dist_vec_init, Re_vec_init);
    plt.plot(init_dist, Re_2, 'go');
    plt.plot(init_dist + dif_dist, Re_3, 'ro');
    plt.plot(init_dist + dif_dist + l_h, Re_4, 'yo');
    plt.plot(init_dist + dif_dist + l_h + dist_vec, Re_vec);
    plt.plot([0, 0.95], [2100, 2100], '-');
    plt.xlabel('Comprimento [m]');
    plt.ylabel('Numero de Reynolds');
    plt.xticks(np.arange(0, 1 ,0.1));
    plt.yticks(np.arange(0, 1.1 * Re_1, 1000))
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.show();
    Re_chart.savefig("/Users/Libotte/Desktop/re_chart.pdf", bbox_inches='tight')
    
    Vel_chart = plt.figure(); #Grafico da velocidade
    
    plt.plot(dist_vec_init, v_vec_init);
    plt.plot(init_dist, v2, 'go');
    plt.plot(init_dist + dif_dist, v3, 'ro');
    plt.plot(init_dist + dif_dist + l_h, v4, 'yo')
    plt.plot(init_dist + dist_vec + dif_dist + l_h, v_vec);
    plt.xlabel('Comprimento [m]');
    plt.ylabel('Velocidade [m/s]');
    plt.xticks(np.arange(0, 1 ,0.1));
    plt.yticks(np.arange(0, v1+0.2, 0.2));
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.show();
    Vel_chart.savefig("/Users/Libotte/Desktop/vel_chart.pdf", bbox_inches='tight')
    
    Pre_chart = plt.figure(); #Grafico da pressao
    
    plt.plot([0, init_dist], [pre_1, pre_2]);
    plt.plot(init_dist, pre_2, 'go');
    plt.plot(init_dist + dif_dist, pre_3, 'ro');
    plt.plot(init_dist + dif_dist + l_h, pre_4, 'yo')
    plt.plot([init_dist + dif_dist + l_h, init_dist + dif_dist + l_h + sec_dist], [pre_4, pre_5]);
    plt.xlabel('Comprimento [m]');
    plt.ylabel('Pressao [Pa]');
    plt.xticks(np.arange(0, 1 ,0.1));
    plt.yticks(np.arange(pre_1 - 13, pre_1+0.5, 1));
    plt.grid(color='k', linestyle='-', linewidth=0.5)
    plt.show();
    Pre_chart.savefig("/Users/Libotte/Desktop/pre_chart.pdf", bbox_inches='tight')
    
    return None


project(x1, x2, etta, angle)





