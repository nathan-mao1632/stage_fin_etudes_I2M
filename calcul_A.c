/* =====================================================
 * INFOS GENERALES
 * =====================================================
 *
 * obj : calculer la matrice de passage A qui et le produit entre les distances au sein des mailles * dmu/dT
 *
 * code qui calcule la trajectoire en 2d de rayons dans un 
 * echantillon avec un maillage (I,J) de réfraction et absorption
 * au sein de chaque maille
 *
 * A prendre en compte : 
 *   - Fonctionne pour inclinaison (alpha) entre 0 et pi
 *   - Choisir le bon champ de mu à comparer pour le calcul du poids (le faire automatiquement mais je l'ai pas fait encore)
 */



#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>


/* =========================================================
   Structures
   ========================================================= */
typedef struct {
    double *data;
    int size;
    int capacity;
} Tableau_dynamique_angles;

typedef struct {
    double *data;
    int size;
    int capacity;
} Tableau_dynamique_points_intersections;

typedef struct {
    int *data;
    int size;
    int capacity;
} Tableau_dynamique_maillage;

typedef struct {
    double *data;
    int size;
    int capacity;
} Tableau_dynamique_distance_mailles;


/* =========================================================
   Init / ajout / affichage / libération
   ========================================================= */
void init_angles(Tableau_dynamique_angles *arr) {
    arr->data = NULL; arr->size = 0; arr->capacity = 0;
}
void init_points_intersections(Tableau_dynamique_points_intersections *arr) {
    arr->data = NULL; arr->size = 0; arr->capacity = 0;
}
void init_maillage(Tableau_dynamique_maillage *arr) {
    arr->data = NULL; arr->size = 0; arr->capacity = 0;
}
void init_distance_mailles(Tableau_dynamique_distance_mailles *arr) {
    arr->data = NULL; arr->size = 0; arr->capacity = 0;
}




static void grow_angles(Tableau_dynamique_angles *arr) {
    if (arr->size + 2 > arr->capacity) {
        int nc = arr->capacity == 0 ? 4 : arr->capacity * 2;
        arr->data = realloc(arr->data, nc * sizeof(double));
        arr->capacity = nc;
    }
}
void ajouter_point_intersection(Tableau_dynamique_points_intersections *arr, double x, double y) {
    if (arr->size + 2 > arr->capacity) {
        int nc = arr->capacity == 0 ? 4 : arr->capacity * 2;
        arr->data = realloc(arr->data, nc * sizeof(double));
        arr->capacity = nc;
    }
    arr->data[arr->size]     = x;
    arr->data[arr->size + 1] = y;
    arr->size += 2;
}

void ajouter_point_maille(Tableau_dynamique_maillage *arr, int x, int y) {
    if (arr->size + 2 > arr->capacity) {
        int nc = arr->capacity == 0 ? 4 : arr->capacity * 2;
        arr->data = realloc(arr->data, nc * sizeof(int));
        arr->capacity = nc;
    }
    arr->data[arr->size]     = x;
    arr->data[arr->size + 1] = y;
    arr->size += 2;
}

void ajouter_distance_maille(Tableau_dynamique_distance_mailles *arr, double x) {
    if (arr->size + 1 > arr->capacity) {
        int nc = arr->capacity == 0 ? 2 : arr->capacity *2;
        arr->data = realloc(arr->data, nc * sizeof(double));
        arr->capacity = nc;
    }
    arr->data[arr->size]     = x;
    arr->size += 1;
}
/*==================================================
 * fonction d'échantillonage de la position d'arrivée sur le pixel (entre Pmin et Pmax)
 * ================================================*/
double randu(double a, double b){
        return a+(b-a)*(double)rand()/RAND_MAX;}

/*========================================================
   Calcul pour la discrétisation (somme des distances au sein des mailles * le coef d'absorption correspondant)
   ======================================================*/

double calcul_discretisation(Tableau_dynamique_maillage *mailles, Tableau_dynamique_distance_mailles *distance, 
		int I, int J, 
		double *Abs_flat, 
		double resultat_discretisation){
	for (int i = 0; i< distance->size; i++) {
        	resultat_discretisation += distance->data[i]
                	*Abs_flat[(mailles->data[i])*J+ mailles->data[i+1]];
	}
	return resultat_discretisation;
}


/* =========================================================
   Fonctions angle (Snell-Descartes)
   ========================================================= */
void calcul_angle_verti_verti(Tableau_dynamique_angles *arr, double alpha,
                               double theta_2, double theta_prime, double n_1, double n_2) {
    grow_angles(arr);
    double t = asin(n_2 * sin(theta_2) / n_1);
    arr->data[arr->size]     = t;
    if ((alpha >= M_PI/2) && (theta_prime<=0)){
        arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha >= M_PI/2) && (theta_prime>0)){
        arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha < M_PI/2) &&(theta_prime<=0)){
        arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha < M_PI/2) &&(theta_prime>0)){
        arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    arr->size += 2;
}
void calcul_angle_hori_verti_interface_haute(Tableau_dynamique_angles *arr, double alpha,
                              double theta_2, double theta_prime, double n_1, double n_2) {
    grow_angles(arr);
    double t = -asin(n_2 * sin(M_PI / 2 - theta_2) / n_1);
    arr->data[arr->size]     = t;
    //pas possible d'avoir H-V haute quand (alpha>=M_PI/2) et (theta_prime<=0)
    if ((alpha>=M_PI/2) && (theta_prime>0)){
       arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime<=0)){
       arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime>0)){
       arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;} /*j'avais un plus avant mais marchait pas pour pi/3*/
    arr->size += 2;
}
void calcul_angle_hori_verti_interface_basse(Tableau_dynamique_angles *arr, double alpha,
                              double theta_2, double theta_prime, double n_1, double n_2) {
    grow_angles(arr);
    double t = asin(n_2 * sin(M_PI / 2 - theta_2) / n_1);
    arr->data[arr->size]     = t;
    if ((alpha>=M_PI/2) && (theta_prime<=0)){
       arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha>=M_PI/2) && (theta_prime>0)){
       arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime<=0)){
       arr->data[arr->size + 1] = -M_PI / 2 + alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime>0)){
       arr->data[arr->size + 1] = M_PI / 2 + alpha - t;}
    arr->size += 2;
}

void calcul_angle_verti_hori_interface_basse(Tableau_dynamique_angles *arr, double alpha,
                              double theta_2, double theta_prime, double n_1, double n_2) {
    grow_angles(arr);
    double t = asin(n_2 * sin(M_PI / 2 - theta_2) / n_1);
    arr->data[arr->size]     = t;
    if ((alpha>=M_PI/2) && (theta_prime<=0)){
       arr->data[arr->size + 1] = -M_PI + alpha - t;}
    if ((alpha>=M_PI/2) && (theta_prime>0)){
       arr->data[arr->size + 1] = -M_PI + alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime<=0)){
       arr->data[arr->size + 1] = -M_PI+ alpha - t;}
    //pas possible d'avoir V-H basse quand (alpha>=M_PI/2) et (theta_prime>=0)
    arr->size += 2;
}
void calcul_angle_verti_hori_interface_haute(Tableau_dynamique_angles *arr, double alpha,
                              double theta_2, double theta_prime, double n_1, double n_2) {
    grow_angles(arr);
    double t = asin(n_2 * sin(M_PI / 2 - theta_2) / n_1);
    arr->data[arr->size]     = t;
    //pas possible d'avoir V-H haute quand (alpha>=M_PI/2) et (theta_prime<=0)
    if ((alpha>=M_PI/2) && (theta_prime>0)){
       arr->data[arr->size + 1] =   alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime<=0)){
       arr->data[arr->size + 1] =  alpha - t;}
    if ((alpha<M_PI/2)&&(theta_prime>0)){
       arr->data[arr->size + 1] =  alpha - t;}
    arr->size += 2;
}
void calcul_angle_hori_hori(Tableau_dynamique_angles *arr, double alpha,
                             double theta_2, double theta_prime, double n_1, double n_2) {
    grow_angles(arr);
    double t = asin(n_2 * sin(theta_2) / n_1);
    arr->data[arr->size]     = t;
    if ((alpha>=M_PI/2) && (theta_prime<=0)){
       arr->data[arr->size + 1] = -M_PI + alpha - t;}
    if ((alpha>=M_PI/2) && (theta_prime>0)){
       arr->data[arr->size + 1] = -M_PI + alpha - t;}
    if ((alpha < M_PI/2) &&(theta_prime<=0)){
        arr->data[arr->size + 1] =  alpha - t;}
    if ((alpha < M_PI/2) &&(theta_prime>0)){
        arr->data[arr->size + 1] = alpha - t;}
    arr->size += 2;
}

void print_array_angles(Tableau_dynamique_angles *arr) {
    printf("[");
    for (int i = 0; i < arr->size; i += 2)
        printf(" (theta_1=%.6f rad, theta_1'=%.6f rad)", arr->data[i], arr->data[i+1]);
    printf(" ]\n\n");
}
void print_array_points(Tableau_dynamique_points_intersections *arr) {
    printf("[");
    for (int i = 0; i < arr->size; i += 2)
        printf(" (%.6e, %.6e)", arr->data[i], arr->data[i+1]);
    printf(" ]\n\n");
}

void print_array_maillage(Tableau_dynamique_maillage *arr) {
    printf("[");
    for (int i = 0; i < arr->size; i += 2)
        printf(" (%d, %d)", arr->data[i], arr->data[i+1]);
    printf(" ]\n\n");
}

void print_array_distance_mailles(Tableau_dynamique_distance_mailles *arr) {
    printf("[");
    double s = 0;
    for (int i = 0; i < arr->size; i += 1){
        printf(" (%.6e)", arr->data[i]);
        s+=arr->data[i];}
    printf(" ]\n\n");
    printf("la somme fait %.6e", s);
}

void free_array_angles(Tableau_dynamique_angles *arr) {
    free(arr->data); arr->data = NULL; arr->size = arr->capacity = 0;
}
void free_array_points_intersections(Tableau_dynamique_points_intersections *arr) {
    free(arr->data); arr->data = NULL; arr->size = arr->capacity = 0;
}
void free_array_maillage(Tableau_dynamique_maillage *arr) {
    free(arr->data); arr->data = NULL; arr->size = arr->capacity = 0;
}
void free_array_distance_mailles(Tableau_dynamique_distance_mailles *arr) {
    free(arr->data); arr->data = NULL; arr->size = arr->capacity = 0;
}



/* =========================================================
   Utilitaires géométriques
   ========================================================= */
static inline double dmin(double a, double b) { return a < b ? a : b; }
static inline double dmax(double a, double b) { return a > b ? a : b; }

/* Tolérance relative pour les tests d'appartenance
   (nécessaire quand alpha=pi/2 et les nœuds ont des erreurs ~1e-18) */
#define EPS_REL 1e-9

/* Retourne 1 si le point (px, py) est sur le segment [P1, P2] */
static int sur_segment(double px, double py,
                        double P1x, double P1y,
                        double P2x, double P2y) {
    double dx = P2x - P1x;
    double dy = P2y - P1y;
    double len2 = dx*dx + dy*dy;
    if (len2 < 1e-30) return 0; /* segment dégénéré */
    /* paramètre t de la projection orthogonale de P sur la droite */
    double t = ((px - P1x)*dx + (py - P1y)*dy) / len2;
    double tol = EPS_REL;
    if (t < -tol || t > 1.0 + tol) return 0;
    /* vérifier que le point est bien sur la droite (pas juste projeté) */
    double qx = P1x + t*dx, qy = P1y + t*dy;
    double dist2 = (px-qx)*(px-qx) + (py-qy)*(py-qy);
    return dist2 < EPS_REL * EPS_REL * len2;
}

/* Intersection rayon / droite en y = y0 + tan_pente*(x - x0)
   Rayon :     y = ym + tan_theta*(xm - x)
   Retourne x de l'intersection. */
static double intersect_V_x(double xm, double ym, double tan_theta,
                           double x0, double y0, double tan_pente) {
    /* Si tan_pente est très grand (interface quasiment verticale x=x0) */
    if (fabs(tan_pente) > 1e14) {
        return x0; /* hypothèse theta=pi/2 donc tan(pi/2)=inf */
    }
    /* Si tan_theta == tan_pente, droites parallèles */
    double denom = tan_pente + tan_theta;
    if (fabs(denom) < 1e-30) return NAN;

    return (ym - y0 + tan_pente * x0 + tan_theta * xm) / denom;
}

static double intersect_H_x(double xm, double ym, double tan_theta,
                           double x0, double y0, double tan_pente) {
    /* Si tan_pente est très grand (interface quasiment verticale x=x0) */
    if (fabs(tan_pente) > 1e14) {
        return x0; /* hypothèse theta=pi/2 donc tan(pi/2)=inf */
    }
    /* Si tan_theta == tan_pente, droites parallèles */
    double denom = -tan_pente + tan_theta;
    if (fabs(denom) < 1e-30) return NAN;

    return (ym - y0 - tan_pente * x0 + tan_theta * xm) / denom;
}


/* =========================================================
   Calcul de trajectoire
   =========================================================
   Le photon remonte depuis le pixel (x_p, y_p) vers la source.
   On traverse le maillage de droite à gauche (on va de J à 0).

   double ***Maillage  : Maillage[i][j] = double[2] avec x=0, y=1
                         indices i = 0..I,  j = 0..J
   double **Refraction : Refraction[i*J+j][0] = indice de la maille (i,j)
   ========================================================= */
void calcul_trajectoire(
        Tableau_dynamique_angles               *angles,
        Tableau_dynamique_points_intersections *pts,
	Tableau_dynamique_maillage *mailles,
	Tableau_dynamique_distance_mailles *distance,
        double ***Maillage,
        int I, int J,
        double alpha,
        double n_air,
        double **Refraction,
	double *Absorption,
        double *Matrice_A_poids_pixel,
        double x_p, double y_p,
	double x_0,
	double theta_echantillonne,
	int *exterieur_maillage,
	double resultat_discretisation)
{
#define N2(i,j)  (Refraction[(i)*J+(j)][0])
#define Mu(i,j)  (Absorption[(i)*J+(j)][0])
#define A(i,j)  (Matrice_A[(i)*(I*J)+(j)][0])

    /* Pentes des deux familles d'interfaces */
    double tan_api  = tan(M_PI - alpha);      /* interfaces verticales   */
    double tan_api2 = tan(fabs(M_PI / 2 - alpha));  /* interfaces horizontales */

    /* --- Point 0 : pixel --- */
    ajouter_point_intersection(pts, x_p, y_p);

    /* --- Point 1 : entrée sur l'interface de droite (colonne J) --- */
    double xref = Maillage[I][J][0], yref = Maillage[I][J][1];
    double x_M1 = intersect_V_x(x_p, y_p, tan(theta_echantillonne), xref, yref, tan_api);
    double y_M1 = y_p + tan(theta_echantillonne) * (x_p - x_M1);

    ajouter_point_intersection(pts, x_M1, y_M1);

    /* Trouver la ligne de la maille d'entrée (dans la colonne J) */
    int i_courant = -1;
    for (int i = 1; i <= I; i++) {
        double ylo = dmin(Maillage[i-1][J][1], Maillage[i][J][1]);
        double yhi = dmax(Maillage[i-1][J][1], Maillage[i][J][1]);
        double tol = EPS_REL * (yhi - ylo);
        if (y_M1 >= ylo - tol && y_M1 <= yhi + tol) { i_courant = i; break; }
    }
    if (i_courant < 0) {
	*exterieur_maillage = 1;
	Matrice_A_poids_pixel[0*(I*J) + 0] = 1.; /*repère pour indiquer que le rayon ne doit pas être pris en compte*/
        /*printf("exterieur_maillage = %d\n", *exterieur_maillage);
        fprintf(stderr, "Erreur : point d'entrée hors maillage (y_M1 = %e)\n", y_M1);
        fprintf(stderr, "  Plage maillage en y : [%e, %e]\n",
                dmin(Maillage[0][J][1], Maillage[I][J][1]),
                dmax(Maillage[0][J][1], Maillage[I][J][1]));*/
        return;
    }
    int j_courant = J;

    /*printf("  Entrée dans maille (%d, %d) au point (%e, %e)\n",
           i_courant, j_courant, x_M1, y_M1);*/

    /* Angles à la première interface (de l'air vers 1ere  maille) */
    double n_maille = N2(i_courant - 1, j_courant - 1);
    double theta_depart = -M_PI/2 - theta_echantillonne + alpha;
/*1ers angles apres la 1ere interface*/
grow_angles(angles);
double t = asin(n_air*sin(theta_depart)/n_maille);
angles->data[angles->size]     =t;
    angles->data[angles->size + 1] = -M_PI/2 + alpha -t; /*pour alpha entre 0 et pi*/
    angles->size += 2;

    //int exterieur_maillage = 0 /* 1 = le rayon est sorti du maillage et ne doit pas être pris en compte*/
    int verticale = 1; /* 1 = dernière interface franchie était verticale */
    int m = 1;         /* indice dans pts du point courant */

    /* =========================================================
       Boucle principale
       ========================================================= */
    while (j_courant > 0) {

        double theta_prime = angles->data[2 * m - 1];
        double tan_theta   = tan(theta_prime);
        double xm = pts->data[2 * m];
        double ym = pts->data[2 * m + 1];
        int i = i_courant;
        int j = j_courant;
        int found = 0;

	

		/*  Interface horizontale du bas */
	    if (alpha>M_PI/2){ /*ce bloc doit être calculé au début que quand alpha>M_PI*/
            if (!found) {
                double xrh = Maillage[i][j-1][0], yrh = Maillage[i][j-1][1];
                double xch= intersect_H_x(xm, ym, tan_theta, xrh, yrh, tan_api2);
                double ych = ym + tan_theta * (xm - xch);

                if (!isnan(xch) &&
                    sur_segment(xch, ych, Maillage[i][j-1][0], Maillage[i][j-1][1], Maillage[i][j][0], Maillage[i][j][1])) {

                    ajouter_point_intersection(pts, xch, ych);
                    ajouter_point_maille(mailles, i_courant, j_courant);
                    double d = sqrt(pow(xch-xm, 2) + pow(ych-ym, 2));
                    ajouter_distance_maille(distance, d);

                    /*Condition de sortie du maillage*/
                    if (i==I){//fprintf(stderr, "Erreur : sortie du maillage (maille (%d, %d))\n", i, j);
                            *exterieur_maillage = 1;
                            //printf("exterieur_maillage = %d\n", *exterieur_maillage);
                            for (int i=0; i<I*J; i++){Matrice_A_poids_pixel[0*(I*J) + i] = 0.;}

                            Matrice_A_poids_pixel[0*(I*J) + 0] = 1.;/*repère pour indiquer que le rayon ne doit pas être pris en compte*/
                    return;}

                    m++;
                    i_courant++;
                    double th_in = angles->data[2*(m-2)];
                    double n1 = N2(i_courant-1, j-1);
                    double n2 = N2(i-1, j-1);
                    if (verticale) calcul_angle_verti_hori_interface_basse(angles, alpha, th_in, theta_prime, n1, n2);
                    else           calcul_angle_hori_hori (angles, alpha, th_in, theta_prime, n1, n2);
                    Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1] = -d*Absorption[(i-1)*J+ j-1];
                    verticale = 0;
                    found = 1;
                }
            }
	}

/*  Interface verticale gauche : colonne j-1  */
            {   
                double xrv = Maillage[i][j-1][0], yrv = Maillage[i][j-1][1];
                double xc;
	        xc= intersect_V_x(xm, ym, tan_theta, xrv, yrv, tan_api);
		
                double yc  = ym + tan_theta * (xm - xc);

		/* interface verticale haute */
                if (!isnan(xc) &&
                    sur_segment(xc, yc, Maillage[i][j-1][0],   Maillage[i][j-1][1], Maillage[i-1][j-1][0], Maillage[i-1][j-1][1])) {
		    ajouter_point_intersection(pts, xc, yc);
		    ajouter_point_maille(mailles, i_courant, j_courant);
		    double d = sqrt(pow(xc-xm, 2) + pow(yc-ym, 2));
		    ajouter_distance_maille(distance, d);
                    m++;

                    /* Condition d'arrêt : on atteint la 1ère interface */
                    if (j -1 == 0) {
                        /*printf("  Trajectoire terminée sur la 1ère interface : (%e, %e)\n",
                               xc, yc);*/
                        j_courant --;
                        double th_in = angles->data[2*(m-2)]; /* theta_1 de la maille précédente */
                        double n1 = n_air;     /* maille d'arrivée */
                        double n2 = N2(i-1, j_courant);       /* maille de départ */
                        if (verticale) calcul_angle_verti_verti(angles, alpha, th_in, theta_prime, n1, n2);
                        else           calcul_angle_hori_verti_interface_haute (angles, alpha, th_in, theta_prime, n1, n2);
                        double theta_new_prime = angles->data[angles->size - 1];
                        double y_source  = yc + tan(theta_new_prime) * (xc - x_0);
                        ajouter_point_intersection(pts, x_0, y_source);
			Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1 ] = -d*Absorption[(i-1)*J+ j-1];
                        m++;
			/*for (int j = 0; j < (I*J); j++) {
                                Matrice_A_poids_pixel[0*(I*J)+j]=0.;}*/
			/*double s=0;
                        for (int i = 0; i < (I*J); i++) {
                                s+=Matrice_A_poids_pixel[0*(I*J)+i];}
                        printf( " Vh  la somme des poids pour un tirage est %.6e\n ",s);
                        printf( "\n");*/
			//for (int j = 0; j < (I*J); j++) {
                          //      Matrice_A_poids_pixel[0*(I*J)+j]=0.;}

			
                        return;
                    }

                    j_courant--;
                    double th_in = angles->data[2*(m-2)]; /* theta_1 de la maille précédente */
                    double n1 = N2(i-1, j_courant-1);     /* maille d'arrivée */
                    double n2 = N2(i-1, j_courant);       /* maille de départ */
                    if (verticale) calcul_angle_verti_verti(angles, alpha, th_in, theta_prime, n1, n2);
                    else           calcul_angle_hori_verti_interface_haute (angles, alpha, th_in, theta_prime, n1, n2);
		    Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1 ] = -d*Absorption[(i-1)*J+ j-1];
		    verticale = 1;
                    found = 1;
                }
	    
		/* interface verticale basse */
                if (!isnan(xc) && i!=I &&
                    sur_segment(xc, yc, Maillage[i+1][j-1][0],   Maillage[i+1][j-1][1], Maillage[i][j-1][0], Maillage[i][j-1][1])) {
                    ajouter_point_intersection(pts, xc, yc);
		    ajouter_point_maille(mailles, i_courant, j_courant);
		    double d = sqrt(pow(xc-xm, 2) + pow(yc-ym, 2));
                    ajouter_distance_maille(distance, d);
                    m++;

                    /* Condition d'arrêt : on atteint la 1ère interface */
                    if (j -1 == 0) {
                        /*printf("  Trajectoire terminée sur la 1ère interface : (%e, %e)\n",
                               xc, yc);*/
                        j_courant --;
                        double th_in = angles->data[2*(m-2)]; /* theta_1 de la maille précédente */
                        double n1 = n_air;     /* maille d'arrivée */
                        double n2 = N2(i-1, j_courant);       /* maille de départ */
                        if (verticale) calcul_angle_verti_verti(angles, alpha, th_in, theta_prime, n1, n2);
                        else           calcul_angle_hori_verti_interface_basse (angles, alpha, th_in, theta_prime, n1, n2);
                        double theta_new_prime = angles->data[angles->size - 1];
                        double y_source  = yc + tan(theta_new_prime) * (xc - x_0);
                        ajouter_point_intersection(pts, x_0, y_source);
                        Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1 ] = -d*Absorption[(i-1)*J+ j-1];
			/*double s=0;
        		for (int i = 0; i < (I*J); i++) {
                		s+=Matrice_A_poids_pixel[0*(I*J)+i];}
        		printf( " Vb  la somme des poids pour un tirage est %.6e\n ",s);
        		printf( "\n");*/
			m++;
			
			return;
		    }

                    j_courant--;
                    double th_in = angles->data[2*(m-2)]; /* theta_1 de la maille précédente */
                    double n1 = N2(i-1, j_courant-1);     /* maille d'arrivée */
                    double n2 = N2(i-1, j_courant);       /* maille de départ */
                    if (verticale) calcul_angle_verti_verti(angles, alpha, th_in, theta_prime, n1, n2);
                    else           calcul_angle_hori_verti_interface_basse (angles, alpha, th_in, theta_prime, n1, n2);
		    Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1] = -d*Absorption[(i-1)*J+ j-1];
		    verticale = 1;
                    found = 1;
                }
            }


            /* Interface horizontale du haut */
            if (!found) {
                double xrh = Maillage[i-1][j-1][0], yrh = Maillage[i-1][j-1][1];
                double xch;
	        if (alpha<=M_PI/2){xch= intersect_V_x(xm, ym, tan_theta, xrh, yrh, tan_api2);}
		else {xch= intersect_H_x(xm, ym, tan_theta, xrh, yrh, tan_api2);}
                double ych = ym + tan_theta * (xm - xch);

                if (!isnan(xch) &&
                    sur_segment(xch, ych, Maillage[i-1][j-1][0], Maillage[i-1][j-1][1], Maillage[i-1][j][0], Maillage[i-1][j][1])) {
		    ajouter_point_intersection(pts, xch, ych);
		    ajouter_point_maille(mailles, i_courant, j_courant);
		    double d = sqrt(pow(xch-xm, 2) + pow(ych-ym, 2));
                    ajouter_distance_maille(distance, d);
		
		    /*Condition de sortie du maillage*/
                    if (i-1==0){//fprintf(stderr, "Erreur : sortie du maillage (maille (%d, %d))\n", i-1, j);
			    i_courant--;
			    *exterieur_maillage = 1;
			    //printf("exterieur_maillage = %d\n", *exterieur_maillage);
			    for (int i=0; i<I*J; i++){Matrice_A_poids_pixel[0*(I*J) + i] = 0.;}
			    Matrice_A_poids_pixel[0*(I*J) + 0] = 1.;/*repère pour indiquer que le rayon ne doit pas être pris en compte*/
                    return;}

                    m++;
                    i_courant--;
                    double th_in = angles->data[2*(m-2)];
                    double n1 = N2(i_courant-1, j-1);
                    double n2 = N2(i-1, j-1);
                    if (verticale) calcul_angle_verti_hori_interface_haute(angles, alpha, th_in, theta_prime, n1, n2);
                    else           calcul_angle_hori_hori (angles, alpha, th_in, theta_prime, n1, n2);
		    Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1 ] = -d*Absorption[(i-1)*J+ j-1];
		    verticale = 0;
                    found = 1;
                }
            }

	        /*  Interface horizontale du bas */
            if (alpha<=M_PI/2){ /*ce bloc doit être calculé à la fin que quand alpha<=M_PI*/
            if (!found) {
                double xrh = Maillage[i][j-1][0], yrh = Maillage[i][j-1][1];
                double xch= intersect_V_x(xm, ym, tan_theta, xrh, yrh, tan_api2);
                double ych = ym + tan_theta * (xm - xch);

                if (!isnan(xch) &&
                    sur_segment(xch, ych, Maillage[i][j-1][0], Maillage[i][j-1][1], Maillage[i][j][0], Maillage[i][j][1])) {

                    ajouter_point_intersection(pts, xch, ych);
                    ajouter_point_maille(mailles, i_courant, j_courant);
                    double d = sqrt(pow(xch-xm, 2) + pow(ych-ym, 2));
                    ajouter_distance_maille(distance, d);

                    /*Condition de sortie du maillage*/
                    if (i==I){//fprintf(stderr, "Erreur : sortie du maillage (maille (%d, %d))\n", i, j);
                            *exterieur_maillage = 1;
                            //printf("exterieur_maillage = %d\n", *exterieur_maillage);
                            for (int i=0; i<I*J; i++){Matrice_A_poids_pixel[0*(I*J) + i] = 0.;}

                            Matrice_A_poids_pixel[0*(I*J) + 0] = 1.;/*repère pour indiquer que le rayon ne doit pas être pris en compte*/
                    return;}

                    m++;
                    i_courant++;
                    double th_in = angles->data[2*(m-2)];
                    double n1 = N2(i_courant-1, j-1);
                    double n2 = N2(i-1, j-1);
                    if (verticale) calcul_angle_verti_hori_interface_basse(angles, alpha, th_in, theta_prime, n1, n2);
                    else           calcul_angle_hori_hori (angles, alpha, th_in, theta_prime, n1, n2);
                    Matrice_A_poids_pixel[0*(I*J) + I*(j-1) + i-1] = -d*Absorption[(i-1)*J+ j-1];
                    verticale = 0;
                    found = 1;
                }
            }
        }




        if (!found) {
            fprintf(stderr,
                    "  Aucune intersection (m=%d, i=%d, j=%d, theta'=%.4f rad)\n",
                    m, i_courant, j_courant, theta_prime);
            return;
        }
        // Sécurité anti-boucle infinie 
        if (m > I * J * 4 + 10) {
            fprintf(stderr, "  Trop d'itérations, arrêt forcé.\n");
            return;
        }
    }


#undef N2
#undef Mu
#undef A
}
/* =========================================================
   Helpers : construire/libérer double*** depuis tableau plat
   ========================================================= */
double*** construire_maillage_ptr(double **src, int rows, int cols) {
    double ***M = malloc(rows * sizeof(double**));
    for (int i = 0; i < rows; i++) {
        M[i] = malloc(cols * sizeof(double*));
        for (int j = 0; j < cols; j++)
            M[i][j] = src[i * cols + j];
    }
    return M;
}
void liberer_maillage_ptr(double ***M, int rows) {
    for (int i = 0; i < rows; i++) free(M[i]);
    free(M);
}

/* =========================================================
   MAIN
   ========================================================= */
int main(void) {

    Tableau_dynamique_angles angles;
    Tableau_dynamique_points_intersections pts;
    Tableau_dynamique_maillage mailles;
    Tableau_dynamique_distance_mailles distance;
    init_angles(&angles);
    init_points_intersections(&pts);
    init_maillage(&mailles);
    init_distance_mailles(&distance);

    FILE *f_infos = fopen("infos.txt", "w");
    FILE *f_chemins = fopen("chemins.dat","w");
    FILE *f_A_moyennes_pixels = fopen("matrice_A_moyennes_pixels.csv","w");
    FILE *f_A_std_pourcent_pixels = fopen("matrice_A_std_pourcent_pixels.csv","w");
    FILE *f_A_poids_pixel = fopen("matrice_A_poids_pixel.csv","w");
    FILE *fichier   = fopen("resultats.dat", "w");
    (void)fichier;

    srand((unsigned)time(NULL));


    /*----- Constantes physiques --- */
    /* double h = 6.62607015e-34, c = 299792458., k = 1.380649e-23; */

    /*Itération Monte Carlo par pixel */
    int N=1;

    /*mesure du temps d'exécution*/
    clock_t start, end;

    /*---- résultats -----*/
    double resultat_discretisation ;

    /*---------champs constants déjà calculés-----------*/
    double champ_constant_200mu_14n_0deg_pisur2 = 0.600000;
    double champ_constant_200mu_14n_0deg_pisur3 = 0.642364;
    double champ_constant_200mu_14n_0deg_pisur4 = 0.695189;

    double champ_constant_200mu_14n_10deg_pisur2 = 0.604669;
    double champ_constant_200mu_14n_10deg_pisur3 = 0.675397;
    double champ_constant_200mu_14n_10deg_pisur4 = 0.739868;

    double champ_choisi_pour_calcul = champ_constant_200mu_14n_10deg_pisur4;

    /* ---- Milieu extérieur ---- */
    double n_air = 1.;

    /* ---- Corps noir ---- */
    /*double theta_1_prime =20. * M_PI / 180; // rad */
    double x_0 = -0.1; /* m */

    /* ---- Capteur ---- */
    double x_p          = 0.1;
    double y_Pmin       = -0.525e-3, y_Pmax = 0.525e-3;/*pour 3 pixels*/
    //double y_Pmin       = -0.2975e-2, y_Pmax = 0.2975e-2;/*pour pi/6*/
    //double y_Pmin       = -0.98e-2, y_Pmax = 0.98e-2;/*pour pi/4*/
    //double y_Pmin       = -2.3e-2, y_Pmax = 2.3e-2; /*pour pi/3*/ 

    double theta_echantillonne =0. * M_PI / 180; /* rad */
    double hauteur_pixel = 0.35e-3; /*m*/
    //double hauteur_pixel = 1.4e-3;  /*m*/

    /*int    nombre_pixels      = 10;
    double hauteur_pixel       = (y_Pmax - y_Pmin) / nombre_pixels;*/

    /* ---- Échantillon ---- */
    double alpha_min    = 1.*M_PI /4 ; /* rad */
    double alpha_max    = 1.*M_PI / 3; /* rad */
    int nombre_rotations    = 0;/*nombre de mesures = nombre_rotations +1 avec la position initiale*/ 
    /* 15 pour avoir une écart de 1° entre pi/4 et pi/3 inclinaison*/


    double mu_0     = 0.2e3;    /* m⁻1 */
    double mu_prime = 0.002;        /* (mK)⁻1 */
    double n_2      = 1.;
    //double longueur = 50e-3;    /* m */
    //double e        = 3e-3;     /* m*/ 


    /* ---- Maillage ---- */
    int I =140, J= 50;

    double longueur = hauteur_pixel*I;
    double e = hauteur_pixel*J;

    /*pour que hauteur_pixel = Maille_verticale*/
    int nombre_pixels = (y_Pmax - y_Pmin)/hauteur_pixel;


    /*Variables utiles dans le code*/
    
    int compteur=0; /*nombre de rayons considérés (ceux qui ne sortent du maillage)*/
    //int *compteur = &cpt;
    int compteur_ligne_A=0; /*numéro de la ligne de A qu'on est en train de remplir (correspond au rayon n - 1)*/
    int ext_maill;
    int *exterieur_maillage = &ext_maill; /*1 = rayon est sorti du maillage donc ne doit pas être pris en compte*/

    /* Tableaux plats de pointeurs */
    double *Abs_flat  = malloc( I * J  * sizeof(double));
    double **Ref_flat  = malloc( I * J * sizeof(double*));
    double *A_poids_pixel_flat  = malloc( (3) * (I*J) * sizeof(double));
    double *A_moyennes_pixels_flat  = malloc( (nombre_pixels*(nombre_rotations+1)) * (I*J) * sizeof(double));
    double *A_std_pixels_flat  = malloc( (nombre_pixels*(nombre_rotations+1)) * (I*J) * sizeof(double));

    /*double mu_par_colonne[J];   // exemple, un indice par colonne
    mu_par_colonne[0] = 200;
    mu_par_colonne[1] = 199;
    mu_par_colonne[2] = 198;
    mu_par_colonne[3] = 197;
    mu_par_colonne[4] = 196;
    mu_par_colonne[5] = 195;
    mu_par_colonne[6] = 194;
    mu_par_colonne[7] = 193;
    mu_par_colonne[8] = 192;
    mu_par_colonne[9] = 191;
    mu_par_colonne[10] =190;*/
    for (int i = 0; i < I; i++) for (int j = 0; j < J; j++) {
        Abs_flat[i*J+j] =-1000.;// mu_par_colonne[j];
    }

    /*double n_par_colonne[J];   // exemple, un indice par colonne
    n_par_colonne[0] = 1.5;
    n_par_colonne[1] = 1.49;
    n_par_colonne[2] = 1.48;
    n_par_colonne[3] = 1.47;
    n_par_colonne[4] = 1.46;
    n_par_colonne[5] = 1.45;
    n_par_colonne[6] = 1.44;
    n_par_colonne[7] = 1.43;
    n_par_colonne[8] = 1.42;
    n_par_colonne[9] = 1.41;
    n_par_colonne[10] = 1.4;*/
    for (int i = 0; i < I; i++) for (int j = 0; j < J; j++) {
        Ref_flat[i*J+j]    = malloc(sizeof(double));
        Ref_flat[i*J+j][0] = n_2;//n_par_colonne[j];
    }



    for (int i = 0; i < (nombre_pixels*(nombre_rotations+1)); i++) for (int j = 0; j < (I*J); j++) {
        A_moyennes_pixels_flat[i*(I*J)+j] = 0.;
    }
    for (int i = 0; i < (nombre_pixels*(nombre_rotations+1)); i++) for (int j = 0; j < (I*J); j++) {
        A_std_pixels_flat[i*(I*J)+j] = 0.;
    }


    /* ---- Affichage des nœuds ---- */
    /*printf(" Nœuds du maillage \n");
    for (int i = 0; i <= I; i++)
        for (int j = 0; j <= J; j++)
            printf("  [%d][%d] = (%10.4e, %10.4e)\n",
                   i, j, M_ptr[i][j][0], M_ptr[i][j][1]);
    printf("\n");*/



/*Boucle sur les angles de rotations de l'échantillons*/

for (int nr = 0; nr<=nombre_rotations; nr++){/*<= car position départ + nombre rotations*/

  double alpha;

  if (nombre_rotations ==0){ alpha = alpha_min;}

   else {alpha = alpha_min + nr*(alpha_max-alpha_min)/(nombre_rotations);}

  double R4[2] = { (-longueur*cos(alpha) - e*sin(alpha)) / 2,
                     ( longueur*sin(alpha) - e*cos(alpha)) / 2 };

  double TrI[2] = { longueur*cos(alpha)/I, -longueur*sin(alpha)/I };
  double TrJ[2] = { e*sin(alpha)/J,  e*cos(alpha)/J};

  double **Mail_flat = malloc((I+1)*(J+1) * sizeof(double*));

  Mail_flat[0*(J+1)+0] = R4;
    for (int i = 1; i <= I; i++) {
        Mail_flat[i*(J+1)+0]    = malloc(2*sizeof(double));
        Mail_flat[i*(J+1)+0][0] = Mail_flat[(i-1)*(J+1)+0][0] + TrI[0];
        Mail_flat[i*(J+1)+0][1] = Mail_flat[(i-1)*(J+1)+0][1] + TrI[1];
    }
    for (int i = 0; i <= I; i++) for (int j = 1; j <= J; j++) {
        Mail_flat[i*(J+1)+j]    = malloc(2*sizeof(double));
        Mail_flat[i*(J+1)+j][0] = Mail_flat[i*(J+1)+(j-1)][0] + TrJ[0];
        Mail_flat[i*(J+1)+j][1] = Mail_flat[i*(J+1)+(j-1)][1] + TrJ[1];
    }

    double ***M_ptr = construire_maillage_ptr(Mail_flat, I+1, J+1);

    /* Écriture pour paraview (tracé pour une inclinaison) */
        fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[0][0][0], M_ptr[0][0][1], 0.);
        fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[0][J][0], M_ptr[0][J][1], 0.);
        fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[I][J][0], M_ptr[I][J][1], 0.);
        fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[I][0][0], M_ptr[I][0][1], 0.);
        fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[0][0][0], M_ptr[0][0][1], 0.);
        fprintf(f_chemins, "\n\n");
        for (int i=1; i<I; i++){
            fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[i][0][0], M_ptr[i][0][1], 0.);
            fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[i][J][0], M_ptr[i][J][1], 0.);
            fprintf(f_chemins, "\n\n");}
        for (int j=1; j<J; j++){
            fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[0][j][0], M_ptr[0][j][1], 0.);
            fprintf(f_chemins, "%.10e %.10e %.10e\n", M_ptr[I][j][0], M_ptr[I][j][1], 0.);
            fprintf(f_chemins, "\n\n");}

    /*  Boucle sur les pixels */
    for (int p = 1; p <= nombre_pixels; p++) {

	/*compter les pixels de bas en haut*/
        /*double y_Pixel_min = y_Pmin + (p-1) * hauteur_pixel;
	double y_Pixel_max = y_Pixel_min + hauteur_pixel;*/

	/*compter les pixels de haut en bas (correspond à ce que fait ASTRA)*/ 
	double y_Pixel_max = y_Pmax - (p-1)* hauteur_pixel;
	double y_Pixel_min = y_Pixel_max - hauteur_pixel;

	for (int i = 0; i < (3); i++) {for (int j = 0; j < (I*J); j++) {
        A_poids_pixel_flat[i*(I*J)+j] = 0.;}
        }
	
	//printf("compteur = %d\n", compteur);

	compteur = 0;  /*nombre de rayons considérés (ceux qui ne sortent du maillage)*/

        //double poids_E=0;     /*poids pour estimer l'espérance*/
        //double poids_V=0;     /*poids pour estimer la variance*/
        //double E, E2, V, S; /* estimateurs de l'espérance, la variance et l'écart-type*/
start = clock();
	for (int n=0; n<N; n++){

        free_array_angles(&angles);
        free_array_points_intersections(&pts);
	free_array_maillage(&mailles);
	free_array_distance_mailles(&distance);
        init_angles(&angles);
        init_points_intersections(&pts);
	init_maillage(&mailles);
	init_distance_mailles(&distance);

	//double w;
	ext_maill = 0; /*1 = rayon est sorti du maillage donc ne doit pas être pris en compte*/

	/*échantillonnage suivant loi uniforme*/
	double y_p = randu(y_Pixel_min, y_Pixel_max);

	/*considère le rayon au milieu du pixel*/
	//double y_p = y_Pixel_min + hauteur_pixel/2;

        //printf(" Pixel %d : (x_p=%e, y_p=%e) \n", p, x_p, y_p);

        calcul_trajectoire(&angles, &pts, &mailles, &distance,
                           M_ptr, I, J, alpha, n_air,
                           Ref_flat, Abs_flat, A_poids_pixel_flat, x_p, y_p, x_0, 
			   theta_echantillonne, &ext_maill,  
			   resultat_discretisation);

	/*double s=0;
	for (int i = 0; i < (I*J); i++) {
		s+=A_poids_pixel_flat[0*(I*J)+i];}
        printf( " la somme des poids pour un tirage est %.6e\n ",s/compteur);
        printf( "\n");*/

	//printf("compteur = %d\n", compteur);
		
	/*printf("Mailles traversées (%d pts) :\n\n", mailles.size / 2);
        print_array_maillage(&mailles);*/ /*
        printf("Points d'intersection (%d pts) :\n\n", pts.size / 2);
        print_array_points(&pts);*/
	/*printf("Distance par mailles (%d distances) :\n\n", distance.size );
        print_array_distance_mailles(&distance);*//*
        printf("Angles (%d paires) :\n\n", angles.size / 2);
        print_array_angles(&angles);*/
	//printf("Le résultat de la discrétisation est %lf\n", calcul_discretisation(&mailles, &distance,I, J, Abs_flat, resultat_discretisation));
       
	compteur_ligne_A++;

	/*Calcul du poids de Monte Carlo*/
/*	if (ext_maill == 1){w=0;}
	else{
	     double a = champ_choisi_pour_calcul;
	     double b = calcul_discretisation(&mailles, &distance,I, J, Abs_flat, resultat_discretisation);
	     w = hauteur_pixel*(a - b)/(1-a);
	     compteur++;}*/


	/*Traçage du rayon*/
        for (int k = 0; k < pts.size; k += 2)
            fprintf(f_chemins, "%.10e %.10e %.10e\n", pts.data[k], pts.data[k+1], 0.);
        //fprintf(fichier, "Le poids est %lf\n", w);
        fprintf(f_chemins, "\n\n");

  /*calcul des estimateurs de la moyenne et écart type par pixel et remplissage de A_moyennes_pixel_flat*/
	
	if ((int)A_poids_pixel_flat[ (0)*(I*J)+0]==1){
		A_poids_pixel_flat[ (0)*(I*J)+0]==0.;}
	else {compteur++;}
	//printf("compteur = %d\n", compteur);
	
	/*double s=0;
        for (int i = 0; i < (I*J); i++) {
                s+=A_poids_pixel_flat[1*(I*J)+i];}
        printf( " la somme des poids pour un tirage est %.6e\n ",s/compteur);
        printf( "\n");*/

	for (int j = 0; j < I*J; j++) {
     			double w = 0.;

             		w = A_poids_pixel_flat[ (0)*(I*J)+j];
			A_poids_pixel_flat[ (0)*(I*J)+j]=0.;

	     /* Incrémentation des estimateurs*/
             		A_poids_pixel_flat[ (1)*(I*J)+j] += w;
             		A_poids_pixel_flat[ (2)*(I*J)+j] += w * w;

	}
	

}/*fin boucle itération Monte Carlo pour 1 pixel*/


/*double s=0;
        for (int i = 0; i < (I*J); i++) {
                s+=A_poids_pixel_flat[1*(I*J)+i];}
        printf( " la somme des poids est %.6e\n ",s);*/
    
end = clock();
//printf("durée du programme : %.6f secondes\n", ((double)end - start) / CLOCKS_PER_SEC);
//printf(" [\n\n");
	for (int j = 0; j < I*J; j++) {

    	  double E = 0.; 
    	  double E2 = 0.;
    	  double V = 0.;
   	  double S = 0.;	
          
	  if (compteur ==0){
	  E = 0.;
          E2 = 0.;
          V = 0.;
          S = 0.;
	  }
	  else {
	  E = A_poids_pixel_flat[ (1)*(I*J)+j]/compteur;
	  E2 = (1.0 /compteur) * A_poids_pixel_flat[ (2)*(I*J)+j];    
	  V = (E2 - E * E)*compteur/(compteur-1);    
	  S = sqrt(V/compteur);
	  }

    	  A_moyennes_pixels_flat[((p)+(nr*nombre_pixels)-1)*(I*J)+j] = E;
    	  A_std_pixels_flat[((p)+(nr*nombre_pixels)-1)*(I*J)+j] = S;

  //        printf( " %.6e; ", A_moyennes_pixels_flat[(p+(nr*nombre_pixels)-1)*(I*J)+j]);

  }

//printf(" ]\n\n");
  
}/*fin boucle pixels*/

liberer_maillage_ptr(M_ptr, I+1);
    for (int i = 1; i <= I; i++) free(Mail_flat[i*(J+1)+0]);
    for (int i = 0; i <= I; i++) for (int j = 1; j <= J; j++) free(Mail_flat[i*(J+1)+j]);
    free(Mail_flat);

}/*fin boucle nombre rotations*/



/*pour avoir l'erreur en pourcentage par rapport à la moyenne*/
for (int i = 0; i < (nombre_pixels*(nombre_rotations+1)); i++) {for (int j = 0; j < (I*J); j++) {
	if (A_moyennes_pixels_flat[i*(I*J)+j] >10e-7){
        A_std_pixels_flat[i*(I*J)+j]=100*A_std_pixels_flat[i*(I*J)+j]/A_moyennes_pixels_flat[i*(I*J)+j];}
}}


 /* printf("Pixel entre %f cm et %f cm\nmu_0 = %f mm-1\nL'angle de départ du faisceau vaut %f radians et alpha vaut %f radians \n\n", y_Pixel_min*pow(10, 2), y_Pixel_max*pow(10, 2), mu_0*pow(10, -3), theta_echantillonne, alpha);

fprintf(fichier, "Pixel entre %f cm et %f cm\nmu_0 = %f mm-1\nL'angle de départ du faisceau vaut %f radians et alpha vaut %f radians \n\n", y_Pixel_min*pow(10, 2), y_Pixel_max*pow(10, 2), mu_0*pow(10, -3), theta_echantillonne, alpha);
//fprintf(fichier, "Les estimateurs ont été calculés avec %d contributions\nL'espérance est estimée à \n%f \navec un écart type de \n%f \n\n", compteur, E, S);
fprintf(fichier, "epaisseur : %f mm, hauteur pixels : %f mm et nombre de pixels : %d \n\n",e*pow(10, 3), hauteur_pixel*pow(10, 3), nombre_pixels);
*/



/* Affichage de tous les rayons du dernier pixel*/
//fprintf(fichier, "\n\n[");
/*for (int i = 0; i < (N); i++) {for (int j = 0; j < (I*J); j++) {
        fprintf(f_A_rayons_pixel, " %.6e ", A_rayons_pixel_flat[i*(I*J)+j]);
    }fprintf(f_A_rayons_pixel, "\n");}*/
//fprintf(fichier, " ]\n\n");

/* Affichage de A avec les estimateurs des espérances sur chaque pixel*/
//fprintf(fichier, "\n\n[");
for (int i = 0; i < (nombre_pixels*(nombre_rotations+1)); i++) {for (int j = 0; j < (I*J); j++) {
	fprintf(f_A_moyennes_pixels, " %.6e; ", A_moyennes_pixels_flat[i*(I*J)+j]);
    }fprintf(f_A_moyennes_pixels, "\n");}
//fprintf(fichier, " ]\n\n");

/* Affichage de A avec les estimateurs des écart-types sur chaque pixel*/
for (int i = 0; i < (nombre_pixels*(nombre_rotations+1)); i++) {for (int j = 0; j < (I*J); j++) {
        fprintf(f_A_std_pourcent_pixels, " %.6e; ", A_std_pixels_flat[i*(I*J)+j]);
    }fprintf(f_A_std_pourcent_pixels, "\n");}

for (int i = 0; i < (nombre_pixels*(nombre_rotations+1)); i++) {
	double s = 0;
	for (int j = 0; j < (I*J); j++) {s+=A_moyennes_pixels_flat[i*(I*J)+j];}
        fprintf(f_infos, " somme ligne : %.6e\n ", s);
    }





fprintf(f_infos, "Construction Matrice_A avec :\n\n inclinaison entre %f et %f degres\n angle %f degres\n indice de refraction : %f\n Maillage : %d par %d\n plage de pixels entre %fcm et %fcm\n nombre de rotations : %d\n pas d'angle : %f\n\n La matrice est calculée avec :\n\n longueur d'echantillon %fmm\n epaisseur %fmm\n une hauteur de pixels de %fmm\n nombre de pixels : %d\n nombre de tirages par pixel N=%d \n ", alpha_min*180/M_PI,alpha_max*180/M_PI, theta_echantillonne*180/M_PI, n_2, I, J, y_Pmin*pow(10, 2), y_Pmax*pow(10, 2), nombre_rotations, (alpha_max-alpha_min)*180/(M_PI*nombre_rotations), longueur*pow(10,3), e*pow(10,3), hauteur_pixel*pow(10,3), nombre_pixels, N);


    /*  Libération */
    free_array_angles(&angles);
    free_array_points_intersections(&pts);
    free_array_maillage(&mailles);
    free_array_distance_mailles(&distance);
    for (int i = 0; i < I; i++) for (int j = 0; j < J; j++) {
        //free(Abs_flat[i*J+j]);
        free(Ref_flat[i*J+j]);
    }
    free(Abs_flat); free(Ref_flat);
    /*for (int i = 0; i < (N*nombre_pixels); i++) for (int j = 0; j < I*J; j++) {
        free(A_flat[i*(I*J)+j]);
    }*/
    free(A_poids_pixel_flat);free(A_moyennes_pixels_flat);
    free(A_std_pixels_flat);
    fclose(f_chemins); fclose(fichier);fclose(f_infos);
    fclose(f_A_poids_pixel);fclose(f_A_moyennes_pixels);
    fclose(f_A_std_pourcent_pixels);


    return 0;
}
