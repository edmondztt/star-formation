#include"define.h"

/*** with conduction, no artificial support ***/

int main()
{
	start = clock();
	int    i = 0, j = 0, k = 0;
	double mytime = 0, TimeContinue = 0, time1;
	int ff = 0, ff2 = 0;
	int countcolumn = 0, countcolumncf = 0;
	double result[15][100], result_star[15][100], Ek_evo[15][100], frag[15][100], coag[15][100];  //result[0] stores the initial distribution; frag[] & coag[] store the separate contributions;

    /***************initialize*************/
    initialize(&mytime, &TimeContinue, &ff, &ff2);
	time1 = (double)mytime/myr;

	/************open file****************/
	char f_txt[25];	
    char f_distr[25];
    char f_mass[25];
    
	char f_star[25];
	char f_binary[25];
	char f_companion[25];
    char cp[25];

	FILE *fp, *fp1, *fpfrag, *fpcoag, *fstar, *fdistr, *fbinary, *fcompanion;
	filename_gen(f_txt, f_binary, f_distr, f_star, f_companion, cp);
	
/*********** below do star formation *************/
	i=99;
	while(m[i]>=mBEnew)
	{
		star_formation(i);
		star_formation_binary(i);
		star_formation_companion(i);
		i--;
	}

///////////*******************************///////////

	finish = clock();

	fstar=fopen(f_star,"w");
	fbinary=fopen(f_binary,"w");
	fdistr=fopen(f_distr,"w");
	fcompanion=fopen(f_companion,"w");
	output_distr(fdistr);
	output_star(fstar);
	output_binary(fbinary);
	output_companion(fcompanion);
	fclose(fdistr);
	fclose(fstar);
	fclose(fbinary);
	fclose(fcompanion);
        return 1;

}

/**************functions***************/
int output_distr(FILE* fdistr)
{
	int i=0;
	for(i=0; i<100; i++)
	{
		fprintf(fdistr,"%lf\t", n[i]);
	}
	return 1;
}

int output_companion(FILE* fcompanion)
{
	int i=0;
	for(i=0; i<100; i++)
	{
		fprintf(fcompanion,"%lf\t%lf\n",n_star_sec[i], n_star_prim[i]);
	}
	return 1;
}

int output_star(FILE* fstar)
{
	int i=0;
	for(i=0; i<100; i++)
	{
		fprintf(fstar,"%lf\t",n_star[i]);
	}
	return 1;
}

int output_binary(FILE* fbinary)
{
	int i=0;
	for(i=0; i<100; i++)
	{
		fprintf(fbinary,"%lf\n",n_binary[i]);
	}
	return 1;
}

int filename_gen(char f_txt[25], char *f_binary, char*f_distr, char*f_star, char* f_companion, char*cp)
{
	int i;
	printf("please input the number of the txt file:\n");
	start2 = clock();
//	scanf("%s", f_txt);
	f_txt[0]='\0';
	strcat(f_txt,"0311-01.txt");
	end2 = clock();
        strcpy(f_distr, f_txt);
	strcpy(f_star, f_txt);
	strcpy(f_binary, f_txt);
	strcpy(f_companion, f_txt);
	strcpy(cp,f_txt);
	i = 0;
	do
	{
		i++;
	}
	while(f_distr[i] != '.');
	do
	{
		f_distr[i] = '\0';
		f_star[i] = '\0';
		f_binary[i] = '\0';
		f_companion[i] = '\0';
		cp[i] = '\0';
		i++;
	}
	while(f_distr[i] != '\0');
	strcat(f_distr, "-distr_r.txt");
	strcat(f_star, "-star.txt");
	strcat(f_binary, "-binary.txt");
	strcat(f_companion, "-companion.txt");
	strcat(cp, "-cp.txt");
	return 1;
}


int finalize(double mytime, FILE *fp, char *cp, int ff, int ff2)
{
	int i;
	double time1;
	time1 = (double)(mytime/myr);
	fprintf(fp, "%lf myr\n", time1);
	fprintf(fp, "cost time: %fs\n", (double)(finish - start + start2 - end2)/CLOCKS_PER_SEC);
	fclose(fp);
	FILE *fcp;
	fcp = fopen(cp,"w");
	fprintf(fcp, "%f\n", mytime);
	fprintf(fcp, "%f\n", star_mass);
	fprintf(fcp, "%f\n", back_mass);
	fprintf(fcp, "%f\n", ghost_mass);
	for(i = 0; i<100; i++)
	{
		fprintf(fcp, "%f\n", n[i]);
	}
	fprintf(fcp, "\n");
	for(i = 0; i<100; i++)
	{
		fprintf(fcp, "%f\n", v_d[i]);
	}
	for(i = 0; i<100; i++)
	{
		fprintf(fcp, "%f\n", n_star[i]);
	}
	fprintf(fcp, "%d\n%d\n", ff, ff2);
	fclose(fcp);
	return 1;
};

int print_paras(FILE *fp, double back_mass, double TimeContinue, double total_mass_core)
{
	float al = alpha, be = beta, q_be = q_beta, sf = SF_rate * myr, g1 = Gauss1, g2 = Gauss2, fif = fifa;
	float BEmass = mBE, rb = Rb/pc, bm = (double)back_mass/Msolar, timestep = Timestep / yr, eva = Evap_rate*100*kyr,drho = Drho, rhoc1=rhob*Drho*1e16;
	if(TimeContinue)
	{
		fprintf(fp,"this continues file %s\ninitial time is: %fmyr\n", name, TimeContinue/myr);
	}
	fprintf(fp, "bin: %d\nalpha: %f\nbeta: %f\nSF_rate: %f /myr\nGaussian sigma1: %f\nGaussian sigema2: %f\nfilling factor: %f\n",bin, al, be, sf, g1, g2, (double)total_mass_core/back_mass/drho);
	fprintf(fp, "Rb: %f pc\nBackground mass: %f Msolar\ntimestep: %f yr\nevaporation rate: %f%% per kyr\nrhoc: %fe-16\nDensity contrast: %f\n", rb, bm, timestep, eva, rhoc1, drho);
	fprintf(fp, "v_m index: %f\n", v_m);
	fprintf(fp, "Bonnor-Ebert mass: %f Msolar\n", BEmass/Msolar);
	fprintf(fp, "dn/dq = q^q_beta, q_beta=%f\n", q_be);
	fprintf(fp, "\ncondensation kappa=%f\n", kappa);
	fprintf(fp, "\nvelocity dispersion=%f\n", vc);
	return 1;
}

int readdata(double *mytime, double *TimeContinue, int *ff, int *ff2)
{
	int i, j;
	char y[25], z[25], name[25];
	FILE *fpdata;
	printf("please input the number of the check point file:\n");
//	scanf("%s", name);
	name[0]='\0';
	strcat(name,"0311-01-cp.txt");
	fpdata = fopen(name, "r");
	if(fpdata==NULL)
	{
		printf("checkpoint file not found!\n");
		return 0;
	}

	fscanf(fpdata, "#mytime = %lf\n", mytime);
	fscanf(fpdata, "#star_mass = %lf\n", &star_mass);
	fscanf(fpdata, "#back_mass = %lf\n", &back_mass);
	fscanf(fpdata, "#ghost_mass = %lf\n", &ghost_mass);
	fscanf(fpdata, "#Nbin N VD NSTAR\n");
	fscanf(fpdata, "#[1] [2] [3] [4]\n");
	for(i = 0; i<100; i++)
	{
		fscanf(fpdata, "%d %lf %lf %lf\n", &j,n+i,v_d+i,n_star+i);
	}

	fclose(fpdata);
	return 1;
}


int initialize(double *mytime, double *TimeContinue, int *ff, int*ff2)
{
	int i, j;
	int init = 10;
	double a = 1.07897;//a=exp((log(100)-log(0.05))/2)=m[i+1]/m[i];
	double A;
	char flag;
	
//constants initialize:
	m[0] = 0.05 * Msolar;
	sum[0] = 0;
	m_insolar[0] = (double)m[0]/Msolar;
	m[100] = 100 * Msolar;
	delta_n[0] = 0;
	delta_n_coag[0] = 0;
	delta_n_frag[0] = 0;
	delta_n_cond[0] = 0;
	delta_n_evap[0] = 0;
	for( i = 1; i < 100; i++)
	{ 
		m[i] = m[i-1] * a;              
		m_insolar[i] = (double)m[i]/Msolar;
		sum[i] = sum[i-1] + m[i-1];
		delta_n[i] = 0;
		delta_n_coag[i] = 0;
		delta_n_frag[i] = 0;
		delta_n_cond[i] = 0;
		delta_n_evap[i] = 0;
		n_star[i] = 0;
		n_binary[i] = 0;
		n_star_sec[i]=0;
		n_star_prim[i]=0;
		for(j = 0; j<4; j++)
		{
			rk[j][i] = 0;
		}
	}
	n_star[0] = 0;
	n_binary[0]=0;
	n_star_sec[0]=0;
	n_star_prim[0]=0;
//interaction:
    printf("Use the default initial settings? y/n:\n");
	flag = 'n';
	if (flag != 'y')
	{
		readdata(mytime, TimeContinue, ff, ff2);
		n_star[0] = 0;/// just for follow 0513-11-cp.txt
		total_mass_core = 0;
		totEk = 0;
		for(i = 0; i<100; i++)
		{
			total_mass_core += n[i]*m[i];
		}
		total_mass = total_mass_core + back_mass + star_mass;
		return 1;
	}

	return 1;
};

int star_formation_hierarchy(int i)
{
	return 1;
}

int star_formation(int i)
{
//	double temp = sf_rate[i] * n[i] * Timestep;
	double temp_star_mass, tempn_star;
//	delta_n[i] -= temp;
//	star_mass += n[i] * alpha * m[i];
	temp_star_mass = m[i] * alpha;
	int j, k;
	j=0;
	k=j;
	while(m[j]<temp_star_mass)
	{
		k=j;
		j++;
	}
	tempn_star = n[i] * alpha * m[i]/m[k];
	n_star[k] += tempn_star;
//  	back_mass += temp * (1 - alpha) * m[i];          /*alpha is star formation efficiency*/
	return 1;
};

int star_formation_companion(int i)
{
	int m1bin=38;  //corresponds to the new mBE~0.83Msolar;
	int j,k;
//	for(j=0; j<100; j++)
//	{ 
//		n_star_sec[j] = 0;
///		n_star_prim[j] = 0;
//	}
	//below do the sf: m1 in the range of (m[0], m[31]<0.5Msolar<m[32]); dn/dm1 = 1/m1;
	double m1=0, m2, n1, n2;
	j=0;
	for(j=0; j<m1bin; j++)
	{
		n1 = n[i] * (pow(m[j+1],lambda_comp+1)-pow(m[j],lambda_comp+1))/(pow(m[m1bin], lambda_comp+1)-pow(m[0],lambda_comp+1));
		n_star_sec[j] += n1;
		n2 = n1;
		m2 = m[i] - m[j];
//		k=m1bin-1;
		k=j;
		while(m[k]<=m2)
			k++;
		k--;
		n_star_prim[k] += n2*m2/m[k];
	}
	return 1;
}

int star_formation_binary(int i)
{
	int j, flag;
	double q, q1, Qb=q_beta+1, tt1, tt2, tt, tmpstar;  //tmpstar for debug;
	double mtmp[2], ntmp;
	double A_q;
	double q0=0.01;
	double dq=0.01;
	tmpstar = 0;
	for(q=q0; q<=0.5-dq; q+= dq)  //pow(0,-0.5) will have some problem??
	{
		mtmp[0] = m[i]*q;
		mtmp[1] = m[i]*(1-q);
		ntmp = n[i];
		A_q=(pow(0.5,Qb)-pow(q0,Qb));
		if(ntmp>0)
		q1 = q+dq;
		tt1 = pow(q1,Qb);
		tt2 = pow(q,Qb);
		tt = (tt1-tt2);
		ntmp *= tt/A_q;
		if(ntmp<0)
			ntmp=-ntmp;
		for(flag=0; flag<2; flag++)
		{
			j=0;
			while(m[j]<=mtmp[flag])
				j++;
			j--;
			if(j>=0)
				n_binary[j] += mtmp[flag]/m[j]*ntmp;
		}
	}
	return 1;
}

