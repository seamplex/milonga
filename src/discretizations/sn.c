/*------------ -------------- -------- --- ----- ---   --       -            -
 *  milonga's discrete ordinates common routines
 *
 *  Copyright (C) 2015--2016 jeremy theler
 *
 *  This file is part of milonga.
 *
 *  milonga is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  milonga is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with wasora.  If not, see <http://www.gnu.org/licenses/>.
 *------------------- ------------  ----    --------  --     -       -         -
 */

#include <gsl/gsl_math.h>

#include "../milonga.h"

// direcciones y pesos para ordenadas discretas
double **Omega;
double *w;

// tabla 4-1 de lewiss (p 162)
// coincide con la tabla 1 de la pag. 208 de stammler-abbate
// TODO: que se pueda elegir para s4 en adelante
#define S2_MU1    1.0/M_SQRT3
#define S2_W1     1.0

#define S4_MU1    0.3500212
#define S4_W1     (1.0/3.0)

#define S6_MU1    0.2666355
#define S6_W1     0.1761263
#define S6_W2     0.1572071

#define S8_MU1    0.2182179
#define S8_W1     0.1209877
#define S8_W2     0.0907407
#define S8_W3     0.0925926

// spot rule
int sn_init_weights(void) {
  
  int m, n, j;
  
  // inicializacion de pesos y direcciones SN
  w = calloc(milonga.directions, sizeof(double));
  Omega = calloc(milonga.directions, sizeof(double *));
  for (m = 0; m < milonga.directions; m++) {
    Omega[m] = calloc(3, sizeof(double));
  }
  
  if (milonga.dimensions == 1) {
    // tabla 3-1 pagina 121 lewiss
    // en una dimension, las direcciones son los ceros de los polinomios de legendre
    // y los pesos son los que hacen que la cuadratura de gauss sea exacta
    switch (milonga.SN) {
      case 2:
        Omega[0][0] = +1.0/M_SQRT3;
        w[0] =         1.0/2.0;

        Omega[1][0] = -Omega[0][0];
        w[1] =         w[0];
      break;
      case 4:
        Omega[0][0] =  sqrt(3.0/7.0 - 2.0/7.0 * sqrt(6.0/5.0));
        w[0] =         0.6521451549/2.0;

        Omega[1][0] =  sqrt(3.0/7.0 + 2.0/7.0 * sqrt(6.0/5.0));
        w[1] =         0.3478548451/2.0;

        Omega[2][0] = -Omega[0][0];
        w[2] = w[0];

        Omega[3][0] = -Omega[1][0];
        w[3] = w[1];
      break;

      case 6:
        Omega[0][0] =  0.2386191860;
        w[0] =         0.4679139346/2.0;

        Omega[1][0] =  0.6612093864;
        w[1] =         0.3607615730/2.0;

        Omega[2][0] =  0.9324695142;
        w[2] =         0.1713244924/2.0;

        Omega[3][0] = -Omega[0][0];
        w[3] =         w[0];

        Omega[4][0] = -Omega[1][0];
        w[4] =         w[1];

        Omega[5][0] = -Omega[2][0];
        w[5] =         w[2];
      break;

      case 8:
        Omega[0][0] =  0.1834346424;
        w[0] =         0.3626837834/2.0;

        Omega[1][0] =  0.5255324099;
        w[1] =         0.3137066459/2.0;

        Omega[2][0] =  0.7966664774;
        w[2] =         0.2223810344/2.0;

        Omega[3][0] =  0.9602898564;
        w[3] =         0.1012285363/2.0;

        Omega[4][0] = -Omega[0][0];
        w[4] =         w[0];

        Omega[5][0] = -Omega[1][0];
        w[5] =         w[1];

        Omega[6][0] = -Omega[2][0];
        w[6] =         w[2];

        Omega[7][0] = -Omega[3][0];
        w[7] =         w[3];
      break;

      default:
        wasora_push_error_message("unsupported N (%d)", milonga.SN);
        return WASORA_RUNTIME_ERROR;
      break;
    }
  } else {
    // tabla 4-1 pagina 162 de lewiss
    // tambien se podria usar el archivo sndir2.dat de fentraco, pero es diferente
    double s4_mu1 = S4_MU1;
    double s4_mu2 = sqrt(S4_MU1*S4_MU1 + (2-6*S4_MU1*S4_MU1) * (2-1) / (4-2));

    double s6_mu1 = S6_MU1;
    double s6_mu2 = sqrt(S6_MU1*S6_MU1 + (2-6*S6_MU1*S6_MU1) * (2-1) / (6-2));
    double s6_mu3 = sqrt(S6_MU1*S6_MU1 + (2-6*S6_MU1*S6_MU1) * (3-1) / (6-2));

    // esto es para 2 y 3 dimensiones    
    int N_octs = (milonga.dimensions == 2) ? 4 : 8;
    int J_octs = milonga.directions / N_octs;

    // ponemos los pesos como todas las permutaciones en el primer cuadrante
    // y despues  rellenamos los otros
    switch (milonga.SN) {
      case 2:
        Omega[0][0] = S2_MU1;
        Omega[0][1] = S2_MU1;
        // hay que analizar por que tenemos que poner esto en cero pero
        // si no lo hacemos falla el scattering anisotropico en 2d
        Omega[0][2] = (milonga.dimensions == 3)?S2_MU1:0;
        w[0] =         S2_W1/(double)(N_octs);
      break;
      case 4:
        Omega[0][0] = s4_mu1;
        Omega[0][1] = s4_mu1;
        Omega[0][2] = (milonga.dimensions == 3)?s4_mu2:0;
        w[0] =         S4_W1/(double)(N_octs);

        Omega[1][0] = s4_mu1;
        Omega[1][1] = s4_mu2;
        Omega[1][2] = (milonga.dimensions == 3)?s4_mu1:0;
        w[1] =         S4_W1/(double)(N_octs);

        Omega[2][0] = s4_mu2;
        Omega[2][1] = s4_mu1;
        Omega[2][2] = (milonga.dimensions == 3)?s4_mu1:0;
        w[2] =         S4_W1/(double)(N_octs);
      break;
      case 6:
        Omega[0][0] = s6_mu1;
        Omega[0][1] = s6_mu1;
        Omega[0][2] = (milonga.dimensions == 3)?s6_mu3:0;
        w[0] =         S6_W1/(double)(N_octs);

        Omega[1][0] = s6_mu1;
        Omega[1][1] = s6_mu2;
        Omega[1][2] = (milonga.dimensions == 3)?s6_mu2:0;
        w[1] =         S6_W2/(double)(N_octs);
        
        Omega[2][0] = s6_mu2;
        Omega[2][1] = s6_mu1;
        Omega[2][2] = (milonga.dimensions == 3)?s6_mu2:0;
        w[2] =         S6_W2/(double)(N_octs);
        
        Omega[3][0] = s6_mu2;
        Omega[3][1] = s6_mu2;
        Omega[3][2] = (milonga.dimensions == 3)?s6_mu1:0;
        w[3] =         S6_W2/(double)(N_octs);

        Omega[4][0] = s6_mu3;
        Omega[4][1] = s6_mu1;
        Omega[4][2] = (milonga.dimensions == 3)?s6_mu1:0;
        w[4] =         S6_W1/(double)(N_octs);
        
        Omega[5][0] = s6_mu1;
        Omega[5][1] = s6_mu3;
        Omega[5][2] = (milonga.dimensions == 3)?s6_mu1:0;
        w[5] =         S6_W1/(double)(N_octs);

        break;
      default: 
        wasora_push_error_message("to be done");
        return WASORA_RUNTIME_ERROR;
      break;
    }
    
    // ya tenemos el primer octante, ahora rellenamos todos los otros
    for (n = 1; n < N_octs; n++) {
      for (j = 0; j < J_octs; j++) {
        Omega[n*J_octs + j][0] = ((n & 1) ? (-1) : (+1)) * Omega[j][0];
        Omega[n*J_octs + j][1] = ((n & 2) ? (-1) : (+1)) * Omega[j][1];
        Omega[n*J_octs + j][2] = ((n & 4) ? (-1) : (+1)) * Omega[j][2];
        w[n*J_octs + j] = w[j];
      }
    }
  }

  return WASORA_RUNTIME_OK;
}