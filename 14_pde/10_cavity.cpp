#include <cstdio>
#include <cmath>
#include <vector>

int main(){
	int nx = 41;
	int ny = 41;
	int nt = 500;
	int nit = 50;
	float dx = 2.0/(nx - 1);
	float dy = 2.0/(ny - 1);
	float dt = 0.01;
	int rho = 1;
	float nu = 0.02;

	std::vector<std::vector<float> > u(ny, std::vector<float>(nx, 0.0));
	std::vector<std::vector<float> > v(ny, std::vector<float>(nx, 0.0));
	std::vector<std::vector<float> > p(ny, std::vector<float>(nx, 0.0));
	std::vector<std::vector<float> > b(ny, std::vector<float>(nx, 0.0));
	
	for(int n = 0; n <= nt; n++){
		for(int j = 1; j <= ny - 2; j++){
			for(int i = 1; i <= nx - 2; i++){
				b[j][i] = rho * (1.0 / dt *\
						((u[j][i+1] - u[j][i-1]) / (2.0 * dx) + (v[j+1][i] - v[j-1][i]) / (2.0 * dy)) -\
						pow((u[j][i+1] - u[j][i-1]) / (2.0 * dx),2) - 2.0 * ((u[j+1][i] - u[j-1][i]) / (2.0 * dy) *\
						(v[j][i+1] - v[j][i-1]) / (2.0 * dx)) - pow((v[j+1][i] - v[j-1][i]) / (2.0 * dy),2.0));
			}	
		}
		
		for(int it = 0; it <= nit; it++){
			std::vector<std::vector<float> > pn(ny, std::vector<float>(nx, 0.0));
			for(int j = 0; j <= nx - 1; j++){
		       		for(int i = 0; i <= ny - 1; i++){
					pn[j][i] = p[j][i];
					}
			}
			for(int j = 1; j <= ny - 2 ; j++){
				for(int i = 1; i <= nx - 2; i++){
					p[j][i] = (pow(dy,2)*(pn[j][i+1] + pn[j][i-1]) +\
							pow(dx,2)*(pn[j+1][i] + pn[j-1][i]) -\
							b[j][i] * pow(dx,2) * pow(dy,2))\
							/ (2 * (pow(dx,2) + pow(dy,2)));
				}
			}
			for(int i = 0; i <= ny - 1; i++){
		       	p[i][ny-1] = p[i][ny-2];
			}
			for(int i = 0; i <= nx - 1; i++){
                p[0][i] = p[1][i];
            }
			for(int i = 0; i <= ny - 1; i++){
                p[i][0] = p[i][1];
            }
			for(int i = 0; i <= ny - 1; i++){
                p[nx-1][i] = 0;
            }
		}
		std::vector<std::vector<float> > un(ny, std::vector<float>(nx, 0.0));
        for(int j = 0; j <= nx - 1; j++){
            for(int i = 0; i <= ny - 1; i++){
                un[j][i] = u[j][i];
            }
		}
		std::vector<std::vector<float> > vn(ny, std::vector<float>(nx, 0.0));
        for(int j = 0; j <= nx - 1; j++){
            for(int i = 0; i <= ny - 1; i++){
                vn[j][i] = v[j][i];
            }
		}
		for(int j = 1; j <= ny - 2; j++){
			for(int i = 1; i <= nx - 2; i++){
				u[j][i] = un[j][i] -un[j][i] * dt / dx *(un[j][i] - un[j][i-1])\
					  - un[j][i] * dt / dy * (un[j][i] - un[j-1][i])\
					  - dt / (2 * rho * dx) * (p[j][i+1] - p[j][i-1])\
					  + nu * dt / pow(dx,2) * (un[j][i+1] - 2 * un[j][i] +un[j][i-1])\
					  + nu * dt / pow(dy,2) * (un[j+1][i] -2 * un[j][i] + un[j-1][i]);

				v[j][i] = vn[j][i] -vn[j][i] * dt / dx *(vn[j][i] - vn[j][i-1])\
                                          - vn[j][i] * dt / dy * (vn[j][i] - vn[j-1][i])\
                                          - dt / (2 * rho * dx) * (p[j+1][i] - p[j-1][i])\
                                          + nu * dt / pow(dx,2) * (vn[j][i+1] - 2 * vn[j][i] +vn[j][i-1])\
                                          + nu * dt / pow(dy,2) * (vn[j+1][i] -2 * vn[j][i] + vn[j-1][i]);
			}
		}

		for(int i = 0; i <= nx - 1; i++){
			u[0][i] = 0;
		}
		for(int i = 0; i <= nx - 1; i++){
            u[i][0] = 0;
        } 
		for(int i = 0; i <= nx - 1; i++){
            u[i][nx-1] = 0;
        }
		for(int i = 0; i <= nx - 1; i++){
            u[nx-1][i] = 1;
        }
		for(int i = 0; i <= nx - 1; i++){
            v[0][i] = 0;
        }
		for(int i = 0; i <= nx - 1; i++){
            v[nx-1][i] = 0;
        }
		for(int i = 0; i <= nx - 1; i++){
            v[i][0] = 0;
        }
		for(int i = 0; i <= nx - 1; i++){
            v[i][nx-1] = 0;
        }
		printf("%f ",v[1][1]);
	}

}






			


