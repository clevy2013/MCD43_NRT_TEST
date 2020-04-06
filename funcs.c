#include <math.h>
#include <assert.h>
#include "def.h"
#include "funcs.h"
#include "util.h"
#include "model.h"
#include "input.h"
#include "log.h"

extern double North_boundary_lat;
extern double South_boundary_lat;

#ifdef DEBUG
int tot_rm = 0;
#endif

extern double bsaLiSparseRker[91];
extern double bsaRossThickker[91];
extern double wsaKer[3];
extern double nadirRef45Ker[3];

static Output *m_output_250m = NULL;
static Output *m_output_500m = NULL;
static Output *m_output_1km = NULL;

void init_proc()
{
#ifdef R250M
	if(m_output_250m == NULL){
		m_output_250m = new_output(0, 0, NCOL_250M);
	}
#else
	if(m_output_500m == NULL){
		m_output_500m = new_output(0, 0, NCOL_500M);
	}

	if(m_output_1km == NULL){
		m_output_1km = new_output(0, 0, NCOL_1KM);
	}
#endif
}

void clean_proc()
{
	if(m_output_250m != NULL){
		del_output(m_output_500m);
	}
	if(m_output_500m != NULL){
		del_output(m_output_500m);
	}
	if(m_output_1km != NULL){
		del_output(m_output_1km);
	}
}

Output *get_proc_output(Row *data_row, int jday)
{
	if(data_row->ncol == NCOL_500M){
		m_output_500m->jday = jday;
		m_output_500m->irow = data_row->irow;
		init_output(m_output_500m, NCOL_500M);
		return m_output_500m;
	}
	else if(data_row->ncol == NCOL_1KM){
		m_output_1km->jday = jday;
		m_output_1km->irow = data_row->irow;
		init_output(m_output_1km, NCOL_1KM);
		return m_output_1km;
	}
	else if(data_row->ncol == NCOL_250M){
		m_output_250m->jday = jday;
		m_output_250m->irow = data_row->irow;
		init_output(m_output_250m, NCOL_250M);
		return m_output_250m;
	}
	else{
		return NULL;
	}
}

Output *ProcessADay(Row *data_row, int jday)
{
	int col;
	double med_szn;

	Output *output = get_proc_output(data_row, jday);
	if(output == NULL || data_row->tot_obs == 0){
		return output;
	}

	// select obs of 16 days for this jday
	select_obs(data_row, jday);

	/* all these functions work on Row */
	
	// calc local solar noon
	calc_lsn(output, jday);

	// calc median szn
	calc_med_szn(data_row);

	// determin season for backup db
	determin_season(data_row);
	
	// do inversion
	Obs *obs_list = NULL;
	Obs *pobs = NULL;
	Obs *pobs2 = NULL;
	Brdf *brdf_bk = NULL;
	int n_sel = 0;
	int b;

	for(col=0; col<data_row->ncol; col++){

#ifdef DEBUG
		if(col != DBG_COL){
			continue;
		}
#endif

		if(data_row->cells[col].nobs < 1){
			continue;
		}

		med_szn = data_row->cells[col].med_szn_terra;
		if(med_szn < 0){
			med_szn = data_row->cells[col].med_szn_aqua;
		}
		
		brdf_bk = get_backup_db(data_row, col);

		obs_list = NULL;
		n_sel = 0;

		// init virtual link and outlier flag
		pobs = data_row->cells[col].obs;
		while(pobs != NULL){
			pobs->next_v = NULL;	
			for(b=0; b<NBAND; b++){
				pobs->is_outlier[b] = 0;
			}

			pobs = pobs->next;
		}

		pobs = data_row->cells[col].obs;
		// collect selected obs
		while(pobs != NULL){
			if(pobs->selected){
				if(obs_list == NULL){
					obs_list = pobs;
					pobs2 = obs_list;
				}
				else{
					pobs2->next_v = pobs;
					pobs2 = pobs;
				}

				n_sel ++;
			}	
			pobs = pobs->next;
		}

		if(n_sel < 1){
			continue;
		}

		output->qa[col].is_snow_brdf = data_row->cells[col].is_snow_brdf;

		/* these next functions work on obs list  based on obs->next_v */

		check_outlier(obs_list);

		prepare_inversion(obs_list);

		adjust_weight(obs_list, jday);

		/* for qa */
		platform_flag(obs_list, output, col);
		
		land_water_mask(obs_list, output, col);
		
		has_high_szn_obs(obs_list, output, col);

		/* do actual inversion */
		do_inversion(obs_list, brdf_bk, output, col, med_szn);	

		/* after inversion */
#ifndef NO_BROADBAND
		broad_brdf(&(output->brdf[col]), data_row->cells[col].is_snow_brdf);
	
		assess_overall_qa(output, col, data_row->cells[col].is_snow_brdf);
#endif

		calc_albedo(&(output->brdf[col]), &(output->albedo[col]), &(output->qa[col]), output->lsn);

#ifdef DEBUG
		printf("RESOLUTION:%d RESULT[ROW=%d, COL=%d]:\n", output->ncol, output->irow, col);
		int band;
		for(band=0; band<NBAND+3; band++){
			printf("  BAND %d BRDF=%.4f %.4f %.4f, WSA=%.4f, BSA=%.4f, NBAR=%.4f, INV_STATUS=%d, IS_FILL=%d, QA=%d, MANQA=%d\n",
						band, output->brdf[col].iso[band], output->brdf[col].vol[band], output->brdf[col].geo[band],
						output->albedo[col].wsa[band], output->albedo[col].bsa[band], output->albedo[col].nbar[band],
						output->brdf[col].inv_status[band], output->brdf[col].is_fill[band],
						band<NBAND?output->qa[col].band_qa[band]:output->qa[col].mandatory_qa[band], output->qa[col].mandatory_qa[band]);
		}
		printf("\n");
#endif
	
		// 2/18/16 sqs
		// to pass to the 1km data_row, for kicking fake snow
		/*if(data_row->ncol == NCOL_500M){
			if(output->is_fill[col]){
				data_row->cells[col].is_snow_brdf = QA_FILL;
			}
			else if(output->qa[col].is_snow_brdf == 1){
				data_row->cells[col].is_snow_brdf = 1;
			}
			else{
				data_row->cells[col].is_snow_brdf = 0;
			}
		}*/

	}	// for col

	do_count(data_row, output);

	return output;
}

int combine_two_500m_rows(Row *row1, Row *row2, Row *row_1km)
{
	int col, i, j;
	int cnt;

	if(row1 == NULL && row2 == NULL){
		return -1;
	}

	Row *rows_500m[2] = {row1, row2};

	Row *prow = row1;
	if(prow == NULL){
		prow = row2;
	}	

	clear_row(row_1km);
	row_1km->irow = prow->irow / 2;

	// check if is two blank rows 
	int tot = 0;
	if(row1 != NULL){
		tot += row1->tot_obs;
	}
	if(row2 != NULL){
		tot += row2->tot_obs;
	}
	if(tot == 0){
		return 0;
	}

	Obs *pobs;
	Obs *pobs2;

	Bkdb *bk_ptr;
	int bk_cnt[3][NBAND];

	int s, b, p;

	for(col=0; col<row_1km->ncol; col++){
#ifdef DEBUG
		if(col != DBG_COL){
			continue;
		}
#endif
		memset(bk_cnt, 0, 3*NBAND*sizeof(int));
		bk_ptr = row_1km->cells[col].brdf_bk;
		float *bk_pars[3][3] = {{bk_ptr->bk_s1->iso, bk_ptr->bk_s1->vol, bk_ptr->bk_s1->geo},
														{bk_ptr->bk_s2->iso, bk_ptr->bk_s2->vol, bk_ptr->bk_s2->geo},
														{bk_ptr->bk_s3->iso, bk_ptr->bk_s3->vol, bk_ptr->bk_s3->geo}};

		// sqs 6/24/15: init 1km backupdb for sumup
		for(s=0; s<3; s++){
			for(p=0; p<3; p++){
				for(b=0; b<NBAND; b++){
					bk_pars[s][p][b] = 0.0;
				}
			}
		}

		pobs = row_1km->cells[col].obs;
		cnt = 0;

		for(i=0; i<2; i++){
			if(rows_500m[i] == NULL){
				continue;
			}

			for(j=0; j<2; j++){
				// 2/18/16 sqs
				// check if any of the 2x2 500m pixels is snow_brdf
				// then give the 1km pixel a initial snow status
				// for kicking out fake snow.
				// new: only 4 500m snow free then change the 1km snow status
				/*if(rows_500m[i]->cells[col*2+j].is_snow_brdf == 0){
					snow_free_cnt ++;
					if(snow_free_cnt > 2){
						row_1km->cells[col].is_snow_brdf = 0;	
					}
				}*/

				pobs2 = rows_500m[i]->cells[col*2+j].obs;
				while(pobs2 != NULL){

					if(pobs == NULL){
						row_1km->cells[col].obs = pobs2;
						pobs = pobs2;
					}
					else{
						pobs->next = pobs2;
						pobs = pobs2;
					}
					cnt ++;
					pobs2 = pobs2->next;
				}

				// average brdf bk
				bk_ptr = rows_500m[i]->cells[col*2+j].brdf_bk;	
				Brdf *seasons[3] = {bk_ptr->bk_s1, bk_ptr->bk_s2, bk_ptr->bk_s3};
				float *bk_pars_500m[3][3] = {{bk_ptr->bk_s1->iso, bk_ptr->bk_s1->vol, bk_ptr->bk_s1->geo},
													 					{bk_ptr->bk_s2->iso, bk_ptr->bk_s2->vol, bk_ptr->bk_s2->geo},
																		{bk_ptr->bk_s3->iso, bk_ptr->bk_s3->vol, bk_ptr->bk_s3->geo}};

				for(s=0; s<3; s++){
					for(b=0; b<NBAND; b++){
						if(!seasons[s]->is_fill[b]){
							for(p=0; p<3; p++){
								bk_pars[s][p][b] += bk_pars_500m[s][p][b];
#ifdef DEBUG
								printf("  ADDING UP BACKUP DB: SEASON=%d, BAND=%d, PAR=%d, BRDF=%f\n",
											 s, b, p, bk_pars_500m[s][p][b]);
#endif
							}
							bk_cnt[s][b] ++;
						}
					}
				}

			}
		}

		row_1km->cells[col].nobs = cnt;
		row_1km->tot_obs += cnt;

		bk_ptr = row_1km->cells[col].brdf_bk;
		Brdf *seasons[3] = {bk_ptr->bk_s1, bk_ptr->bk_s2, bk_ptr->bk_s3};
		for(s=0; s<3; s++){
			for(b=0; b<NBAND; b++){
				if(bk_cnt[s][b] > 0){
					for(p=0; p<3; p++){
						bk_pars[s][p][b] /= bk_cnt[s][b];
					}
					seasons[s]->is_fill[b] = 0;
				}
				else{
					seasons[s]->is_fill[b] = 1;
				}
			}
		}

	}	// col loop

	// break the link
	if(row1 != NULL){
		for(col=0; col<row1->ncol; col++){
			row1->cells[col].obs = NULL;
		}
	}
	if(row2 != NULL){
		for(col=0; col<row2->ncol; col++){
			row2->cells[col].obs = NULL;
		}
	}

	return 0;
}

void calc_med_szn(Row *data_row)
{
	int col;
	Obs *pobs;
	int i;
	int j;
	float szns_terra[10000];
	float szns_aqua[10000];
	int cnt_terra;
	int cnt_aqua;
	float *szns;
	int *cnt;

	for(col=0; col<data_row->ncol; col++){

		if(data_row->cells[col].nobs < 1){
			data_row->cells[col].med_szn_terra = DBL_FILL;
			data_row->cells[col].med_szn_aqua = DBL_FILL;
			continue;
		}

		cnt_terra = 0;
		cnt_aqua = 0;
		pobs = data_row->cells[col].obs;
		while(pobs != NULL){	
			if(!pobs->selected){
				pobs = pobs->next;
				continue;
			}

			if(pobs->platform == MODIS_TERRA){
				szns = szns_terra;
				cnt = &cnt_terra;
			}
			else{
				szns = szns_aqua;
				cnt = &cnt_aqua;
			}

			// insert into the ordered-list
			for(i=0; i<*cnt; i++){
				if(pobs->szn < szns[i]){
					for(j=*cnt-1; j>i; j--){
						szns[j] = szns[j-1];
					}		
					break;
				}
			}	
			szns[i] = pobs->szn;
			*cnt = *cnt + 1;
				
			pobs = pobs->next;
		}

		data_row->cells[col].med_szn_terra = szns_terra[cnt_terra/2];
		data_row->cells[col].med_szn_aqua = szns_aqua[cnt_aqua/2];
	}
}

Brdf * get_backup_db(Row *data_row, int icol)
{
	Brdf *pbrdf = NULL;

	int season = data_row->cells[icol].season;
	switch(season){
		case SEASON_ACTIVE:
			pbrdf = data_row->cells[icol].brdf_bk->bk_s1;
			break;
		case SEASON_INACTIVE:
			pbrdf = data_row->cells[icol].brdf_bk->bk_s2;
			break;
		case SEASON_SNOW:
			pbrdf = data_row->cells[icol].brdf_bk->bk_s3;
			break;
	}	

	return pbrdf;
}

void determin_season(Row *data_row)
{
	int i;
	double ndvi;
	int cnt;
	Obs *pobs;

	for(i=0; i<data_row->ncol; i++){
		if(data_row->cells[i].nobs < 1){
			continue;
		}

		if(data_row->cells[i].is_snow_brdf == 1){
			data_row->cells[i].season = SEASON_SNOW;
		}
		else{
			ndvi = 0.0;
			cnt = 0;
			pobs = data_row->cells[i].obs;
			while(pobs != NULL){
				if(pobs->selected){
					ndvi += (pobs->ref[1] - pobs->ref[0])/(pobs->ref[1] + pobs->ref[0]);
					cnt ++;
				}
				pobs = pobs->next;
			}

			if(cnt>0){
				ndvi /= (double)cnt;
			}

			if(ndvi < NDVI_ACTIVE){
				data_row->cells[i].season = SEASON_INACTIVE;
			}
			else{
				data_row->cells[i].season = SEASON_ACTIVE;
			}
		}
	}
}

void check_outlier(Obs *obs_list)
{
	Obs *pobs;
	double mean[NBAND];
	double diff[NBAND];
	double sd[NBAND];
	double cv[NBAND];
	double max_sd[NBAND];
	int cnt[NBAND];
	double tmp;
	int b;

	for(b=0; b<NBAND; b++){
		mean[b] = 0.0;
		diff[b] = 0.0;
		sd[b] = 0.0;
		cv[b] = 0.0;
		cnt[b] = 0;
	}

	// calc mean	
	pobs = obs_list;
	while(pobs != NULL){
		for(b=0; b<NBAND; b++){
			if(!pobs->is_fill[b]){
				mean[b] += pobs->ref[b];
				cnt[b] ++;
			}
		}

		pobs = pobs->next_v;
	}

	for(b=0; b<NBAND; b++){
		if(cnt[b] > 0){
			mean[b] /= (double)cnt[b];
		}
		else{
			mean[b] = 0.0;
		}
	}

	// calc sd, cv
	pobs = obs_list;
	while(pobs != NULL){
		for(b=0; b<NBAND; b++){
			if(!pobs->is_fill[b]){
				diff[b] += (mean[b] - pobs->ref[b]) * (mean[b] - pobs->ref[b]);
			}
		}

		pobs = pobs->next_v;
	}
	for(b=0; b<NBAND; b++){
		sd[b] = sqrt (diff[b] / (double) (cnt[b] - 1));
		cv[b] = sd[b] * 100.0 / mean[b];
	}
	
	// outlier 
	for(b=0; b<NBAND; b++){
		if(cnt[b] <= 2){
			continue;
		}

		if ((cv[b] + 0.0001) > MAX_CV){
			max_sd[b] = -1.0;
			pobs = obs_list;
			while(pobs != NULL){
				if(pobs->is_fill[b]){
					pobs = pobs->next_v;
					continue;
				}

				tmp = fabs(pobs->ref[b] - mean[b]) / sd[b];
				if(tmp > max_sd[b]){
					max_sd[b] = tmp;
				}
				
				pobs = pobs->next_v;
			}

			pobs = obs_list;
			while(pobs != NULL){
				if(pobs->is_fill[b]){
					pobs = pobs->next_v;
					continue;
				}

				tmp = fabs(pobs->ref[b] - mean[b]) / sd[b];
				if (tmp >= (max_sd[b] - 0.0001)){
					pobs->is_outlier[b] = 1;
#ifdef DEBUG
					printf("Outlier: band=%d, ref=%f, mean=%f, cv=%f, sd=%f\n", b, pobs->ref[b], mean[b], cv[b], sd[b]);
#endif
				}
				
				pobs = pobs->next_v;
			}
		}
	}

}

void adjust_weight(Obs *obs_list, int jday_interest)
{
	int b;
	Obs *pobs;
	pobs = obs_list;
	
#ifdef NRT
	/*double Lap_w[16] = { 1, 0.986666667, 0.973333333, 0.96, 0.946666667, 0.933333333, 0.92, 0.906666667,
											 0.893333333, 0.88, 0.866666667, 0.853333333, 0.84, 0.826666667, 0.813333333, 0.8 };*\
	double Lap_w[16] = { 1, 0.99, 0.98, 0.97, 0.96, 0.95, 0.94, 0.93, 
				0.92, 0.91, 0.90, 0.89, 0.88, 0.87, 0.85, 0.85 };
#else
	double Lap_w[9] = { 1, 0.9167, 0.8404, 0.7704, 0.7062, 0.6474, 0.5935, 0.5441, 0.4987 };	/*Laplace distribution */
#endif

	int dist;

	while(pobs != NULL){
		dist = abs(pobs->jday - jday_interest);

		for(b=0; b<NBAND; b++){
			if(!pobs->is_fill[b]){
				pobs->weight_v[b] = pobs->weight[b] * Lap_w[dist];
			}
		}

		pobs = pobs->next_v;
	}
}

/*pobs->weight_v[b] = 0.4 * pobs->weight[b] + 0.3 * Lap_w[dist] + 0.3 * pobs->obscov;*\
/*pobs->weight_v[b] = 0.5 * pobs->weight[b] + 0.5 * Lap_w[dist];*\

/* calc kernels for selected obs 
 * weight may change after each selection
 * so this function should be called 
 * after select_obs() */
void prepare_inversion(Obs *obs_list)
{
	Obs *pobs;
	pobs = obs_list;
	while(pobs != NULL){
		CalculateKernels(pobs->kerval, pobs->vzn*D2R, pobs->szn*D2R, pobs->ran*D2R);

		pobs = pobs->next_v;
	}
}

void do_inversion(Obs *obs_list, Brdf *brdf_bk, Output *output, int col, double med_szn)
{
	int band;
	int nobs;
	Obs *pobs;
	int n_fill_band = 0;

	int jday = output->jday;

	for(band=0; band<NBAND; band++){
		nobs = 0;
		pobs = obs_list;
		while(pobs != NULL){
			if(!pobs->is_fill[band] && !pobs->is_outlier[band]){
				nobs ++;
			}
			pobs = pobs->next_v;
		}

		if(nobs >= INVNOBSMIN){
			output->brdf[col].inv_status[band] = FULL_INVERSION;
		}
		else if(nobs > 1){
			output->brdf[col].inv_status[band] = MAGN_INVERSION;
		}
		else{
			output->brdf[col].inv_status[band] = FILL_INVERSION;
		}

#ifdef DEBUG
		printf("DO_INVERSION: NOBS=%d, GO %d.\n", nobs, output->brdf[col].inv_status[band]);
#endif

		if(FULL_INVERSION == output->brdf[col].inv_status[band]){

			inversion(obs_list, output, col, band, nobs, jday, med_szn);

			// check rmse, wod
			if(output->rmse[col].band[band]     > BAND_POOR_RMSE ||
				 output->wd_ref45[col].band[band] > BAND_NBAR_POOR_WD ||
				 output->wd_wsa[col].band[band]   > BAND_WSA_POOR_WD){

				output->brdf[col].inv_status[band] = MAGN_INVERSION;
#ifdef DEBUG
				printf("FLIP TO MAGNITUDE INVERSION DUE TO RMSE AND WOD.\n");
#endif
			}

			// check zero pars
			if(EQ(output->brdf[col].vol[band], 0) &&
				 EQ(output->brdf[col].geo[band], 0)){
				
				output->brdf[col].inv_status[band] = MAGN_INVERSION;
#ifdef DEBUG
				printf("FLIP TO MAGNITUDE INVERSION DUE TO ZERO PARAMETERS.\n");
#endif
			}
		}

		if(MAGN_INVERSION == output->brdf[col].inv_status[band]){
			mag_inversion(obs_list, brdf_bk, output, col, band, nobs, jday);
		}

		if(FILL_INVERSION == output->brdf[col].inv_status[band]){
			output->brdf[col].is_fill[band] = 1;	
		}

		// after inersions
		valid_brdf(&(output->brdf[col]), band);

		if(FILL_INVERSION == output->brdf[col].inv_status[band]){
			n_fill_band ++;
		}

		// count valid obs
		count_valid_obs(obs_list, output, col, band);
	
		// band qa
		assess_band_qa(output, col, band, nobs);
		
#ifdef DEBUG
		printf("BAND %d FINAL BRDF: %f, %f, %f, BAND_QA=%d, IS_FILL=%d\n\n", band, output->brdf[col].iso[band],
										output->brdf[col].vol[band], output->brdf[col].geo[band], output->qa[col].band_qa[band],
										output->brdf[col].is_fill[band]);
#endif
	}

	// mandatory qa
	if(n_fill_band == NBAND){
		output->is_fill[col] = 1;
	}
	else{
		output->is_fill[col] = 0;
	}

}

void inversion(Obs *obs_list, Output *output, int col, int band, int nobs, int jday, double med_szn)
{
	int i, j, k;
	double K[3][nobs];
	double R[3] = {0.0, 0.0, 0.0};
	double Matrix[3][3];
	double tmp_matrix[3][3];	// for constrain
	double wt[nobs];
	double pars[3];
	
	int index[3]={0,0,0};
	double d = 0.0;
	double vv[3] = {0.0,0.0,0.0};

	Obs *pobs = obs_list;

	i = 0;
	while(pobs != NULL){
		if(pobs->is_fill[band] || pobs->is_outlier[band]){
			pobs = pobs->next_v;
			continue;	
		}

		wt[i] = pobs->weight_v[band];

		for(j=0; j<3; j++){
			K[j][i] = pobs->kerval[j];
			R[j] += pobs->kerval[j] * pobs->ref[band] * wt[i];
		}

		pobs = pobs->next_v;
		i++;
	}
		
	for(j=0; j<3; j++){
		pars[j] = R[j];
	}

	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			k = 0;
			Matrix[i][j] = 0.0;
			pobs = obs_list;
			while(pobs != NULL){
				if(pobs->is_fill[band] || pobs->is_outlier[band]){
					pobs = pobs->next_v;
					continue;	
				}

				Matrix[i][j] += K[i][k] * K[j][k] * wt[k];
				pobs = pobs->next_v;
				k++;
			}

			// for constrain
			tmp_matrix[i][j] = Matrix[i][j];	
		}
	}
	
	//
	if(0 != ludcmp(Matrix, 3, index, &d, vv)){
		sprintf(msg, "Error! LUDCMP failed.\n");
		log_err();
		exit(-1);
	}
		
	lubksb(Matrix, 3, index, pars);
	
	//
	output->brdf[col].iso[band] = (float)pars[0];
	output->brdf[col].vol[band] = (float)pars[1];
	output->brdf[col].geo[band] = (float)pars[2];
	output->brdf[col].is_fill[band] = 0;

	// constrain
	constrainMeEasy(output, obs_list, col, band, tmp_matrix, med_szn);
	
	pars[0] = output->brdf[col].iso[band];
	pars[1] = output->brdf[col].vol[band];
	pars[2] = output->brdf[col].geo[band];
	
	// calc rmse & WoD
	output->rmse[col].band[band] = 0.0;
	output->wd_ref45[col].band[band] = 0.0;
	output->wd_bsa45[col].band[band] = 0.0;
	output->wd_wsa[col].band[band] = 0.0;

	double modref;
	double diff;

	i=0;
	pobs = obs_list;
	while(pobs != NULL){
		if(pobs->is_fill[band] || pobs->is_outlier[band]){
			pobs = pobs->next_v;
			continue;	
		}

		modref = 0.0;
		for(k=0; k<3; k++){
			modref += pars[k] * K[k][i];
		}
		diff = (modref - pobs->ref[band]) * (modref - pobs->ref[band]);
		output->rmse[col].band[band] += diff * pobs->weight_v[band];

		pobs = pobs->next_v;
		i++;
	}

	double tmp_M[3][3];
	for(i=0; i<3; i++){
		for(j=0; j<3; j++){
			tmp_M[i][j] = 0.0;
			if(i==j){
				tmp_M[i][j] = 1.0;
			}
		}
		lubksb(Matrix, 3, index, tmp_M[i]);
	}

	double wd_ref_buf;
	double wd_bsa_buf;
	double wd_wsa_buf;
	double bsa_ker;

	for(i=0; i<3; i++){
		wd_ref_buf = 0.0;	
		wd_bsa_buf = 0.0;	
		wd_wsa_buf = 0.0;	
		for(j=0; j<3; j++){
			if(j==0){
				bsa_ker = 1.0;
			}
			else if(j==1){
				bsa_ker = bsaRossThickker[45];
			}
			else if(j==2){
				bsa_ker = bsaLiSparseRker[45];
			}

			wd_ref_buf += tmp_M[j][i] * nadirRef45Ker[j];
			wd_bsa_buf += tmp_M[j][i] * bsa_ker;
			wd_wsa_buf += tmp_M[j][i] * wsaKer[j];
		}

		if(i==0){
			bsa_ker = 1.0;
		}
		else if(i==1){
			bsa_ker = bsaRossThickker[45];
		}
		else if(i==2){
			bsa_ker = bsaLiSparseRker[45];
		}
		output->wd_ref45[col].band[band] += wd_ref_buf * nadirRef45Ker[i];
		output->wd_bsa45[col].band[band] += wd_bsa_buf * bsa_ker;
		output->wd_wsa[col].band[band] += wd_wsa_buf * wsaKer[i];
	}

	if(output->wd_ref45[col].band[band] > 0.0){
		output->wd_ref45[col].band[band] = sqrt(output->wd_ref45[col].band[band]);
	}
	else{
		output->wd_ref45[col].band[band] = 0.0;
	}

	if(output->wd_bsa45[col].band[band] > 0.0){
		output->wd_bsa45[col].band[band] = sqrt(output->wd_bsa45[col].band[band]);
	}
	else{
		output->wd_bsa45[col].band[band] = 0.0;
	}

	if(output->wd_wsa[col].band[band] > 0.0){
		output->wd_wsa[col].band[band] = sqrt(output->wd_wsa[col].band[band]);
	}
	else{
		output->wd_wsa[col].band[band] = 0.0;
	}

	if(output->rmse[col].band[band] > 0.0){
		output->rmse[col].band[band] = sqrt(output->rmse[col].band[band]/(double)(nobs-3));
	}
	else{
		output->rmse[col].band[band] = 0.0;
	}

#ifdef DEBUG
	printf("FULL_INV RSLT: row=%d col=%d, band=%d, nobs=%d, jday=%d, pars=%f, %f, %f\n", 
									output->irow, col, band, nobs, jday, pars[0], pars[1], pars[2]);
	printf("               rmse=%f, wd_ref45=%f, wd_bsa45=%f, wd_wsa=%f\n", 
									output->rmse[col].band[band], output->wd_ref45[col].band[band],
									output->wd_bsa45[col].band[band], output->wd_wsa[col].band[band]);
#endif

}

void mag_inversion(Obs *obs_list, Brdf *brdf_bk, Output *output, int col, int band, int nobs, int jday)
{
	double sum1, sum2;
	double modref;
	int i;
	Obs *pobs;
	double pars[3];

	if(brdf_bk == NULL){
		output->brdf[col].is_fill[band] = 1;
		output->brdf[col].inv_status[band] = FILL_INVERSION;
		sprintf(msg, "Error! Backup brdf is null! col=%d, band=%d\n", col, band);
		log_err();
		return;
	}

#ifdef DEBUG
	printf("ENTER MAG_INV: BRDF_DB=%f, %f, %f, IS_FILL=%d\n",
									brdf_bk->iso[band], brdf_bk->vol[band], brdf_bk->geo[band],
									brdf_bk->is_fill[band]);
#endif

	if(brdf_bk->is_fill[band]){
		output->brdf[col].inv_status[band] = FILL_INVERSION;
		return;
	}

	pars[0] = brdf_bk->iso[band];
	pars[1] = brdf_bk->vol[band];
	pars[2] = brdf_bk->geo[band];

	sum1 = 0.0;
	sum2 = 0.0;

	pobs = obs_list;
	while(pobs != NULL){
		if(pobs->is_fill[band] || pobs->is_outlier[band]){
			pobs = pobs->next_v;
			continue;	
		}

		modref = 0.0;
		for(i=0; i<3; i++){
			modref += pobs->kerval[i] * pars[i];
		}

#ifdef DEBUG
		printf("  modref=%f, ref=%f\n", modref, pobs->ref[band]);
#endif

		sum1 += pobs->ref[band] * modref * pobs->weight_v[band];
		sum2 += modref * modref * pobs->weight_v[band];
	
		pobs = pobs->next_v;
	}

	output->brdf[col].iso[band] = brdf_bk->iso[band] * sum1 / sum2;	
	output->brdf[col].vol[band] = brdf_bk->vol[band] * sum1 / sum2;	
	output->brdf[col].geo[band] = brdf_bk->geo[band] * sum1 / sum2;	
	output->brdf[col].is_fill[band] = 0;

#ifdef DEBUG
	printf("MAG_INV RELT: band=%d, iso=%f, vol=%f, geo=%f, nobs=%d\n", band, output->brdf[col].iso[band],
									output->brdf[col].vol[band], output->brdf[col].geo[band], nobs);
	printf("              Backup: band=%d, iso=%f, vol=%f, geo=%f, sum1=%f, sum2=%f\n", band, brdf_bk->iso[band],
									brdf_bk->vol[band], brdf_bk->geo[band], sum1, sum2);
#endif
}

int ProvideObsRow(Row *data_row, Input_files *input, int irow)
{
	// read
	if(0 != ReadObsRow(data_row, input, irow)){
		return -1;
	}

	// prepare
	// this function has been moved into the ReadObsOneRow 
	//PrepareObsRow(data_row);
	
	return 0;
}

void PrepareObsRow(Row *row)
{
	int col;
	Obs *pobs;
	Obs *pobs_pre;
	int flag;

	for(col=0; col<row->ncol; col++){

		pobs_pre = NULL;
		pobs = row->cells[col].obs;

		while(pobs != NULL){
			flag = 0;

			// SZN
			if(pobs->szn >= SZN_POOR){
				flag = 1;			
			}
			// QAs
			else if(0 != check_qa(pobs)){
				flag = 1;			
			}

			if(flag){
				// delete this obs
				if(pobs_pre == NULL){
					row->cells[col].obs = pobs->next;
					del_obs(pobs);
					pobs = row->cells[col].obs;
				}
				else{
					pobs_pre->next = pobs->next;
					del_obs(pobs);
					pobs = pobs_pre->next;
				}

				row->cells[col].nobs --;
				row->tot_obs --;
			}
			else{
				pobs_pre = pobs;
				pobs = pobs->next;
			}

		}
	}
}

// return 0: passed
// return non-zero: drop
int check_qa(Obs *obs)
{
	int MeritValue = 0;
	int qa;//, qa2;
	int l;

	uint16 state_1km = obs->state_1km;
	uint32 QC_500m = obs->QC_500m;

	/* clouds: accept only clear = 00 */
	qa = state_1km & 3;
	if(qa != 0){	//0x11
		return -1;
	}

	/* Feng - (02/01) accept shadows but with punishment */
	qa = state_1km & 4;
	if(qa != 0){	//0x100
		MeritValue += 1;
	}

	/* deep ocean */
	qa = state_1km >> 3;
	qa = qa & 7;						//0x111
	obs->land_water_mask = qa;
	if(qa == SHELFORDEEPOCEAN || qa == FILLVALUELANDWATER){
		return -2;
	}

	/* reject aerosol high = 11, punish climatology = 00 */
	qa = state_1km >> 6;
	qa = qa & 3;						//0x11
	if(qa == 3){						//0x11
		return -3;
	}
	else if(qa == 0){				//0x00
		MeritValue += 1;
	}

	/* reject cirrus high = 11, punish average = 10 */
	qa = state_1km >> 8;
	qa = qa & 3;						//0x11
	if(qa == 3){						//0x11
		return -4;
	}
	else if(qa == 2){				//0x10
		MeritValue += 1;
	}
	
	/* use internal cloud mask as well - Feng (10/04) */
	qa = state_1km >> 10;
	qa = qa & 1;						//0x1
	if(qa == 1){						//0x1
		return -5;
	}
			
	/* snow detected */
	qa = state_1km >> 12;
	qa = qa & 1;						//0x1
	//qa2 = state_1km >> 15;
	//qa2 = qa2 & 1;					//0x1
	//if(qa == 1 || qa2 == 1){
	/* use MOD35 snow mask (bit 12)only - Zhuosen (Jan 2016) */
	if(qa == 1){
		obs->is_snow = 1;
	}
	else{
		obs->is_snow = 0;
	}

	/* BRDF correction */
	// 13 && 14 ??
	
	if (MeritValue > MAX_VALID_MERIT)
		return -6;

	// check individual bands
	for (l = 0; l < NBAND; l++)
	{
		if (MeritValue == 1)
			obs->weight[l] = 0.75;
		if (MeritValue == 2)
			obs->weight[l] = 0.5;
		if (MeritValue > 2 && MeritValue <= MAX_VALID_MERIT)
			obs->weight[l] = 0.25;
	}

	/// band qa
	
	/* cloud effact all bands */
	qa = QC_500m & 3;				//0x11
	if(qa == 2){						//0x10
		return -7;
	}
			
	for (l = 0; l < NBAND; l++)
	{
		MeritValue = 0;

		/* and punish less than ideal quality: 01 */
		qa = QC_500m & 3;				//0x11
		if(qa == 1){						//0x01
			MeritValue += 1;
		}
	
#ifdef R250M
		qa = QC_500m >> (4 + l*4);
#else
		qa = QC_500m >> (2 + l*4);
#endif

		qa = qa & 15;						//0x1111
		if(qa != 0 && qa != 12 && qa != 8){
			MeritValue += 100;
		}
		if(qa == 8){
			MeritValue += 100;
		}
				
		if (MeritValue == 1)
			obs->weight[l] = MIN (0.75, obs->weight[l]);
		if (MeritValue == 2)
			obs->weight[l] = MIN (0.5, obs->weight[l]);
		if (MeritValue > 2 && MeritValue <= MAX_VALID_MERIT)
			obs->weight[l] = MIN (0.25, obs->weight[l]);
		if (MeritValue > MAX_VALID_MERIT){
			obs->is_fill[l] = 1;
			obs->ref[l] = DBL_FILL;
		}
	}

	return 0;
}

/* follow the mod10 algorithem to determin if snow */
/* http://modis-snow-ice.gsfc.nasa.gov/uploads/pap_ESC04GAR.pdf */
/* return 0, 1, 255 (unknown) */
int check_if_snow(Obs *obs)
{
	if(!obs->is_fill[1] && obs->ref[1] <= 0.1){
		return 0;
	}

	if(!obs->is_fill[3] && obs->ref[3] <= 0.1){
		return 0;
	}

	if(obs->is_fill[3] || obs->is_fill[5]){
		return 255;
	}

	double ndsi = (obs->ref[3] - obs->ref[5]) / (obs->ref[3] + obs->ref[5]);
	if(ndsi < 0.1){
		return 0;
	}
	
	if(ndsi >= 0.4){
		return 1;
	}
	
	if(obs->is_fill[0] || obs->is_fill[1]){
		return 255;
	}

	double ndvi = (obs->ref[1] - obs->ref[0]) / (obs->ref[1] + obs->ref[0]);

	if(ndvi + 0.5*ndsi < 0.3){
		return 0;
	}

	if(ndvi > 2.5*ndsi){
		return 0;
	}

	return 1;
}

/* calc snow fraction from ndsi */
/* http://landweb.nascom.nasa.gov/QA_WWW/forPage/snow_user_guide_C6_revised_draft_4.pdf */
uint8 calc_fsc(uint8 ndsi)
{
	if(ndsi > 100){
		return 0;
	}

	double _ndsi = (double)ndsi /100.0;
	int fsc = (int)((-0.01 + (1.45*_ndsi))*100.0);
	if(fsc > 100){
		fsc = 100;
	}
	if(fsc < 0){
		fsc = 0;
	}

	return (uint8)fsc;
}

/* calc obs snow fraction. should be called after check_qa() */
void calc_obs_snow_frac(Obs *obs)
{
	if(!obs->is_snow){
		obs->snow_frac = 0;
		return;
	}

#ifdef R250M
	obs->snow_frac = 100;
	return;	// assume full snow? now use mod10a1
#endif

	// assume snow free by default
	obs->snow_frac = 0;

	if(obs->is_fill[3] || obs->is_fill[5]){
		return;
	}

	uint8 ndsi = (uint8)((obs->ref[3] - obs->ref[5])/(obs->ref[3] + obs->ref[5]) * 100.0);

	if(ndsi > 100){
		ndsi = 100;
	}

	obs->snow_frac = calc_fsc(ndsi);
}


int is_unknown_snow(Obs *obs)
{
	// if band 4 or band 6 is fill, ndsi is unknown
	if(obs->is_snow){
		if(obs->snow_frac == 0 || obs->snow_frac > 100){
			return 1;
		}
	}

	return 0;
}

int is_full_snow(Obs *obs)
{
	if(obs->is_snow && obs->snow_frac >= 95 &&
		 obs->snow_frac <= 100){	// fsc > 95%
		return 1;
	}

	return 0;
}

int is_part_snow(Obs *obs)
{
	if(obs->is_snow && !is_full_snow(obs) && !is_unknown_snow(obs)){
		return 1;
	}

	return 0;
}

/* select obs for this jday (date and snow status)
 * selected obs are flagged for inversion
 */
void select_obs(Row *data_row, int jday)
{
	int i;
	Obs *pobs;
	int nobs;

	// select 16 days based on the jday interest
	for(i=0; i<data_row->ncol; i++){
		nobs = 0;
		pobs = data_row->cells[i].obs;
		while(pobs != NULL){

#ifdef NRT
			if(pobs->jday >= jday-15 && pobs->jday <= jday){
#else
			if(pobs->jday >= jday-8 && pobs->jday <= jday+7){
#endif
				pobs->selected = 1;
				nobs ++;
#ifdef DEBUG
				if(DBG_COL == i){
					int year, doy;
					julday2doy(pobs->jday, &year, &doy);
					printf("SELECTED. DOY=%d, IS_SNOW=%d, SNOW_FRAC=%d RED=%f.\n",
								 doy, pobs->is_snow, pobs->snow_frac, pobs->ref[0]);
				}
#endif
			}
			else{
				pobs->selected = 0;
#ifdef DEBUG
				if(DBG_COL == i){
					int year, doy;
					julday2doy(pobs->jday, &year, &doy);
					printf("UNSELECTED. DOY=%d, IS_SNOW=%d, SNOW_FRAC=%d RED=%f.\n",
							 doy, pobs->is_snow, pobs->snow_frac, pobs->ref[0]);
				}
#endif
			}

			pobs = pobs->next;
		}

#ifdef DEBUG
		if(DBG_COL == i){
			printf("16DAY WINDOW SELECTED OBS: %d\n", nobs);
		}
#endif

		data_row->cells[i].nobs = nobs;
	}	
				
	// snow status
	int dist;
	int min_dist;
	//int num_snow, num_snow_free;
	int num_snow_free, num_full_snow, num_part_snow;
	int num_snow_unkw;
	//double ndsi;
	//double mean_ndsi;
	//int ndsi_cnt;

	for(i=0; i<data_row->ncol; i++){
		// 		get closest obs to day of interest
		min_dist = 16;
		pobs = data_row->cells[i].obs;
		while(pobs != NULL){
			if(!pobs->selected){
				pobs = pobs->next;
				continue;
			}

			dist = abs(pobs->jday - jday);
			if(dist < min_dist){
				min_dist = dist;
				if(min_dist == 0){
					break;
				}
			}
			pobs = pobs->next;
		}

		//   count snow or non-snow, full or partial
		//num_snow = 0;
		num_snow_free = 0;
		num_full_snow = 0;
		num_part_snow = 0;
		num_snow_unkw = 0;

		pobs = data_row->cells[i].obs;
		while(pobs != NULL){
			if(!pobs->selected){
				pobs = pobs->next;
				continue;
			}
			
			dist = abs(pobs->jday - jday);
			if(dist != min_dist){
				pobs = pobs->next;
				continue;
			}

#ifdef DEBUG
		if(DBG_COL == i){
			printf("CHECKING SNOW STA. JDAY=%d, INTEREST=%d, IS_SNOW=%d, SNOW_FRAC=%d\n",
				pobs->jday, jday, pobs->is_snow, pobs->snow_frac);
		}
#endif

			if(!pobs->is_snow){
				num_snow_free ++;
			}
			else if(is_full_snow(pobs)){
				num_full_snow ++;
			}
			else if(is_part_snow(pobs)){
				num_part_snow ++;
			}
			else{
				num_snow_unkw ++;
			}

			pobs = pobs->next;
		}

		int snow_flag;
		if(num_snow_free >= num_full_snow && num_snow_free >= num_part_snow &&
			 num_snow_free >= num_snow_unkw){
			snow_flag = 0;	// snow free
		}
		else if(num_full_snow >= num_snow_free && num_full_snow >= num_part_snow &&
						num_full_snow >= num_snow_unkw){
			snow_flag = 1;	// full snow
		}
		else if(num_part_snow >= num_snow_free && num_part_snow >= num_full_snow &&
						num_part_snow >= num_snow_unkw){
			snow_flag = 2;	// part snow
		}
		else{
			snow_flag = 3;	// unknown snow
		}

		// 2/18/16 sqs
		// if none of the 2x2 500m snow pixels is snow
		// make the 1km pixel snow free
		/*if(snow_flag != 0 && data_row->ncol == NCOL_1KM){
			if(data_row->cells[i].is_snow_brdf == 0){
				snow_flag = 0;
				//debug
				printf("XXXXX Converting 1km pixel (%d, %d) to snow free.\n",
							 data_row->irow, i);
				fflush(stdout);
			}
		}*/
			
		if(snow_flag  == 0){
			data_row->cells[i].is_snow_brdf = 0;
		}
		else{
			data_row->cells[i].is_snow_brdf = 1;
		}

#ifdef DEBUG
		if(DBG_COL == i){
			printf("NUM_NON_SNOW=%d, FRAC_SNOW=%d, FULL_SNOW=%d, SNOW_FLAG=%d, SNOW_BRDF=%d\n",
						 num_snow_free, num_part_snow, num_full_snow, snow_flag, data_row->cells[i].is_snow_brdf);
		}
#endif

		//    remove 
		pobs = data_row->cells[i].obs;
		while(pobs != NULL){
			if(!pobs->selected){
				pobs = pobs->next;
				continue;
			}

			if(snow_flag == 0){
				if(pobs->is_snow){
#ifdef DEBUG
					if(DBG_COL == i){
						printf("REMOVE SNOW OBS: REF=%f\n", pobs->ref[1]);
					}
#endif
					pobs->selected = 0;
					pobs = pobs->next;
					continue;
				}
			}
			else if(snow_flag == 1){
				if(!is_full_snow(pobs) &&
					 !is_unknown_snow(pobs)){
#ifdef DEBUG
					if(DBG_COL == i){
						printf("REMOVE NONE FULL SNOW OBS: REF=%f\n", pobs->ref[1]);
					}
#endif
					pobs->selected = 0;
					pobs = pobs->next;
					continue;
				}
			}
			else if(snow_flag == 2){
				if(!is_part_snow(pobs) &&
					 !is_unknown_snow(pobs)){
#ifdef DEBUG
					if(DBG_COL == i){
						printf("REMOVE NONE PARTIAL SNOW OBS: REF=%f\n", pobs->ref[1]);
					}
#endif
					pobs->selected = 0;
					pobs = pobs->next;
					continue;
				}
			}
			else if(snow_flag == 3){	// keep all snow
				if(!pobs->is_snow){
#ifdef DEBUG
					if(DBG_COL == i){
						printf("REMOVE SNOW FREE OBS: REF=%f\n", pobs->ref[1]);
					}
#endif
					pobs->selected = 0;
					pobs = pobs->next;
					continue;
				}
			}

			pobs = pobs->next;
		}

		// statistics based on nir band
		double mean = 0.0;
		double std = 0.0;
		int cnt = 0;
		pobs = data_row->cells[i].obs;
		while(pobs != NULL){
			if(!pobs->selected){
				pobs = pobs->next;
				continue;
			}

			mean += pobs->ref[1];
			cnt ++;
			pobs = pobs->next;
		}
		
		mean /= (double)cnt;
	
		pobs = data_row->cells[i].obs;
		while(pobs != NULL){
			if(!pobs->selected){
				pobs = pobs->next;
				continue;
			}

			std += (mean - pobs->ref[1]) * (mean - pobs->ref[1]);
			pobs = pobs->next;
		}

		std /= (double)cnt;
		std = sqrt(std);

		if (std > STD_POOR && mean > 0.1 && cnt > 1){
			pobs = data_row->cells[i].obs;
			while(pobs != NULL){
				if(!pobs->selected){
					pobs = pobs->next;
					continue;
				}

				if(pobs->land_water_mask == LAND && pobs->jday != jday){
					pobs->selected = 0;
#ifdef DEBUG
					if(DBG_COL == i){
						printf("REMOVE POOR STD OBS: REF=%f\n", pobs->ref[1]);
					}
#endif
				}
				pobs = pobs->next;
			}
		}

		// the mean ref of the day of interest
		double mean_doi = 0.0;
		if(min_dist == 0){
			cnt = 0;
			pobs = data_row->cells[i].obs;
			while(pobs != NULL){
				if(!pobs->selected){
					pobs = pobs->next;
					continue;
				}

				dist = abs(pobs->jday - jday);
				if(dist != min_dist){
					pobs = pobs->next;
					continue;
				}

				mean_doi += pobs->ref[1];
				cnt ++;
				pobs = pobs->next;
			}
		
			if(cnt > 0){
				mean_doi /= (double)cnt;
				pobs = data_row->cells[i].obs;
				while(pobs != NULL){
					if(!pobs->selected){
						pobs = pobs->next;
						continue;
					}
					if(pobs->ref[1] > 3.0 * mean_doi || pobs->ref[1] < 0.33 * mean_doi){
						if(pobs->ref[1] > 0.02 && pobs->land_water_mask == LAND){
							pobs->selected = 0;
#ifdef DEBUG
							if(DBG_COL == i){
								printf("REMOVE POOR STD2 OBS: REF=%f, MEAN=%f, cnt=%d, MIN_DIST=%d\n", 
												pobs->ref[1], mean_doi, cnt, min_dist);
							}
#endif
						}
					}
					pobs = pobs->next;
				}
			}
		}

	}	// for col
	
}

/*void trim_obs(Row *data_row, int jday)
{
	int i;
	Obs *pobs;

	// select 16 days based on the jday interest
	for(i=0; i<data_row->ncol; i++){
		pobs = data_row->cells[i].obs;
		while(pobs != NULL){
			if(pobs->selected == 1){

			}

			pobs = pobs->next;
		}

	}	
}*/

int ludcmp (double a[3][3], int n, int *indx, double *d, double *vv)
{
   int i;
   int imax=0;
   int j;
   int k;
   double big=0;
   double dum = 0;
   double sum =0;
   double temp = 0;

   *d = 1.0;
   for (i = 0; i < n; i++) {
      big = 0.0;
      for (j = 0; j < n; j++)
				if ((temp = fabs (a[i][j])) > big)
					big = temp;
      if (EQ (big, 0.0)) {
	 			return (-1);
      }
      vv[i] = 1.0 / big;
   }
   for (j = 0; j < n; j++) {
      for (i = 0; i < j; i++) {
	 sum = a[i][j];
	 for (k = 0; k < i; k++)
	    sum -= a[i][k] * a[k][j];
	 a[i][j] = sum;
      }
      big = 0.0;
      for (i = j; i < n; i++) {
	 sum = a[i][j];
	 for (k = 0; k < j; k++)
	    sum -= a[i][k] * a[k][j];
	 a[i][j] = sum;
	 if ((dum = vv[i] * fabs (sum)) >= big) {
	    big = dum;
	    imax = i;
	 }
      }
      if (j != imax) {
	 for (k = 0; k < n; k++) {
	    dum = a[imax][k];
	    a[imax][k] = a[j][k];
	    a[j][k] = dum;
	 }
	 *d = -(*d);
	 vv[imax] = vv[j];
      }
      indx[j] = imax;

  /*if (EQ (a[j][j], 0.0))*/
      if (fabs(a[j][j]) < TINY){
	a[j][j] = TINY;}

      if (j != n - 1) {
	 dum = 1.0 / (a[j][j]);
	 for (i = j + 1; i < n; i++)
	    a[i][j] *= dum;
      }
   }

   return (0);
}

int lubksb (double a[3][3], int n, int *indx, double b[])
{

   int i, ii, ip, j;
   double sum = 0;

   ii = -1;

   for (i = 0; i < n; i++) {
      ip = indx[i];
      sum = b[ip];
      b[ip] = b[i];
      if (ii != -1)  
	 for (j = ii; j <= i - 1; j++)
	    sum -= a[i][j] * b[j];
      else if (sum)
	 ii = i;
      b[i] = sum;
   }

   for (i = n - 1; i >= 0; i--) {
      sum = b[i];
      for (j = i + 1; j < n; j++)
	 sum -= a[i][j] * b[j];
     
       if (fabs(a[i][i]) < TINY ){ /* WZS 1/14 */
	 a[i][i] = TINY;}

        b[i] = sum / a[i][i]; 
     
   }

   return (0);
}

/* CONSTRAIN */

int paramsTest(Brdf *brdf, int band)
{
	float *pars[3];
	pars[0] = brdf->iso;
	pars[1] = brdf->vol;
	pars[2] = brdf->geo;

	int i;
	int ret = 0;

	for(i=0; i<3; i++){
		if(pars[i][band] < -EPS2 || pars[i][band] > 10.0){
			ret = -1;
			break;
		}
		else if(pars[i][band] < 0.0){
			pars[i][band] = 0.0;
		}
	}

#ifdef DEBUG
	printf("PARAME TEST: %s%s%s BRDF=%f %f %f\n", 
									ret<0?KRED:KGRN, ret<0?"FAIL":"PASS", KNRM, pars[0][band], pars[1][band], pars[2][band]);
#endif
	return ret;
}

int albTest(Brdf *brdf, int band)
{
	float *pars[3];
	pars[0] = brdf->iso;
	pars[1] = brdf->vol;
	pars[2] = brdf->geo;

	int i;
	double wsa = 0.0;
	double bsa1 = 0.0;
	double bsa2 = 0.0;

	double bsaKer1[3] = {1.0, bsaRossThickker[0], bsaLiSparseRker[0]};
	double bsaKer2[3] = {1.0, bsaRossThickker[75], bsaLiSparseRker[75]};

	for(i=0; i<3; i++){
		wsa += pars[i][band] * wsaKer[i];
		bsa1 += pars[i][band] * bsaKer1[i];
		bsa2 += pars[i][band] * bsaKer2[i];
	}

	int ret = 0;

	if(wsa<0.0001||wsa>0.9999){
		ret = -1;
	}
	else if(bsa1<0.0001||bsa1>0.9999){
		ret = -2;
	}
	else if(bsa2<0.0001||bsa2>0.9999){
		ret = -3;
	}

#ifdef DEBUG
	printf("ALBEDO TEST: %s%s%s BRDF=%f %f %f, WSA=%f, BSA_LOW=%f BSA_HIGH=%f\n", 
									ret<0?KRED:KGRN, ret<0?"FAIL":"PASS", KNRM, pars[0][band], pars[1][band], pars[2][band], wsa, bsa1, bsa2);
#endif

	return ret;
}

int nadReflTest(Brdf *brdf, int band, double szn)
{
	double kerval[3];	
	CalculateKernels(kerval, 0.0, szn*D2R, 0.0);

	float *pars[3];
	pars[0] = brdf->iso;
	pars[1] = brdf->vol;
	pars[2] = brdf->geo;
	
	int i;
	double nbar = 0.0;
	for(i=0; i<3; i++){
		nbar += kerval[i] * pars[i][band];	
	}

	int ret = 0;
	if(nbar<0.0001){
		ret = -1;
	}

#ifdef DEBUG
	printf("NADIRF TEST: %s%s%s BRDF=%f %f %f, BRF=%f\n", 
									ret<0?KRED:KGRN, ret<0?"FAIL":"PASS", KNRM, pars[0][band], pars[1][band], pars[2][band], nbar);
#endif

	return ret;
}

int testFunction(Brdf *brdf, int band, double szn)
{
	if(0 != paramsTest(brdf, band)){
		return -1;
	}

	if(0 != albTest(brdf, band)){
		return -2;
	}

	if(0 != nadReflTest(brdf, band, szn)){
		return -3;
	}

	return 0;
}

double multiIt (double *x, double *y, int n)
{
	int i;
	double tot = 0.;

	for (i = 0; i < n; i++)
		tot += x[i] * y[i];
	return (tot);
}

void lewisMe (double *_p, double *c, double mat[3][3], int n, int d, double *p)
{
	double _c[4], _pc, _cc;
	double tmp1, tmp2;
	int i, ok;

	/* check that the constraint vector are not 0 */
	ok = 0;
	for (i = 0; i < n; i++)
	{
		if (c[i] != 0)
		{
			ok = 1;
			break;
		}
	}

	if (ok == 1)
	{
		/* calculate C` = M^-1 . C */
		for (i = 0; i < n; i++)
			_c[i] = multiIt (mat[i], c, n);

		/* calculate P` . C */
		_pc = multiIt (_p, c, n);

		/* calculate C` . C */
		_cc = multiIt (_c, c, n);

		for (i = 0; i < n; i++)
		{
			tmp1 = _pc - d;
			if (_cc == 0.)
				tmp2 = 0.;
			else
				tmp2 = (tmp1 / _cc) * _c[i];
			p[i] = _p[i] - tmp2;
		}
	}
	else
	{
		for (i = 0; i < n; i++)
			p[i] = _p[i];
	}
}

void invert_matrix(double mat[3][3], double mat_v[3][3])
{
	int n = 3;
	int index[n];
	double d = 0.0;
	double vv[n];

	int i, j;

	for(i=0; i<n; i++){
		index[i] = 0;
		vv[i] = 0.0;
	}

	ludcmp(mat, n, index, &d, vv);

	for(i=0; i<n; i++){
		for(j=0; j<n; j++){
			mat_v[i][j] = 0.0;
			if(i==j){
				mat_v[i][j] = 1.0;
			}
		}

		lubksb(mat, n, index, mat_v[i]);
	}
}

double calc_rmse(Obs *obs_list, double pars[3], int band)
{
	int i;
	Obs *pobs;
	double modref;
	double diff;
	double rmse = 0.0;

	pobs = obs_list;
	while(pobs != NULL){
		if(pobs->is_fill[band] || pobs->is_outlier[band]){
			pobs = pobs->next_v;
			continue;	
		}

		modref = 0.0;
		for(i=0; i<3; i++){
			modref += pars[i] * pobs->kerval[i];
		}
		diff = (modref - pobs->ref[band]) * (modref - pobs->ref[band]);
		rmse += diff * pobs->weight_v[band];

		pobs = pobs->next_v;
	}

	return rmse;
}

void constrainMeEasy(Output *output, Obs *obs_list, int col, int band, double matrix[3][3], double med_szn)
{
	Brdf *brdf = &(output->brdf[col]);

#ifdef DEBUG
	printf("CONSTRAIN: BRDF=%f, %f, %f, band=%d, median_szn=%.2f\n", brdf->iso[band],
									brdf->vol[band], brdf->geo[band], band, med_szn);
#endif

	// first check
	if(0 == testFunction(brdf, band, med_szn)){
		return;
	}	

	//
	double mat_inv[3][3];

	invert_matrix(matrix, mat_inv);

	double pars[3] = {brdf->iso[band], brdf->vol[band], brdf->geo[band]};
	

	double constrain_kers[6][3] = {{1.0, bsaRossThickker[0],  bsaLiSparseRker[0] }, {1.0, bsaRossThickker[0] , bsaLiSparseRker[0]},	
																 {1.0, bsaRossThickker[75], bsaLiSparseRker[75]}, {1.0, bsaRossThickker[75], bsaLiSparseRker[75]},
																 {1.0, wsaKer[1],           wsaKer[2]          }, {1.0, wsaKer[1],           wsaKer[2]          }};
	int i;
	double p[6][3];
	Brdf brdf_p;
	double rmse;
	double min_rmse = 9999.9;
	int d[6] = {0,1,0,1,0,1};

	int ivalid = -1;

	for(i=0; i<6; i++){
		lewisMe(pars,constrain_kers[i], mat_inv, 3, d[i], p[i]);

		brdf_p.iso[band] = p[i][0];
		brdf_p.vol[band] = p[i][1];
		brdf_p.geo[band] = p[i][2];

		if(0 == testFunction(&brdf_p, band, med_szn)){
			rmse = calc_rmse(obs_list, p[i], band);
			if(rmse < min_rmse){
				ivalid = i;
				min_rmse = rmse;
			}		
		}
	}

	if(ivalid >= 0){
		brdf->iso[band] = p[ivalid][0];	
		brdf->vol[band] = p[ivalid][1];	
		brdf->geo[band] = p[ivalid][2];	
		return;
	}

	output->qa[col].quality_term[band] = 1;

	//
	double conVector[3][3] = {{0.0, 0.0, 1.0}, {0.0, 1.0, 0.0}, {1.0, 0.0, 0.0}};
	min_rmse = 9999.9;
	ivalid = -1;
	for(i=0; i<3; i++){
		lewisMe(pars, conVector[i], mat_inv, 3, 0, p[i]);

		brdf_p.iso[band] = p[i][0];
		brdf_p.vol[band] = p[i][1];
		brdf_p.geo[band] = p[i][2];

		if(0 == testFunction(&brdf_p, band, med_szn)){
			rmse = calc_rmse(obs_list, p[i], band);
			if(rmse < min_rmse){
				ivalid = i;
				min_rmse = rmse;
			}		
		}
	}
	
	if(ivalid >= 0){
		brdf->iso[band] = p[ivalid][0];	
		brdf->vol[band] = p[ivalid][1];	
		brdf->geo[band] = p[ivalid][2];	
		return;
	}

	//
	brdf->iso[band] = 0.0;
	brdf->vol[band] = 0.0;
	brdf->geo[band] = 0.0;

	Obs *pobs = obs_list;
	int cnt = 0;

	while(pobs != NULL){
		if(pobs->is_fill[band] || pobs->is_outlier[band]){
			pobs = pobs->next_v;
			continue;	
		}

		brdf->iso[band] += pobs->ref[band];
		cnt ++;

		pobs = pobs->next_v;
	}

	brdf->iso[band] /= (float)cnt;
	return;
}

void valid_brdf(Brdf *brdf, int band)
{
	int i;	
	float *pars[3] = {brdf->iso, brdf->vol, brdf->geo};

	for(i=0; i<3; i++){
		if(!valid_par(&(pars[i][band]), 0)){
			brdf->is_fill[band] = 1;
			brdf->inv_status[band] = FILL_INVERSION;
			//pars[i][band] = DBL_FILL;
		}
	}
}

void count_valid_obs(Obs *obs_list, Output *output, int col, int band)
{
	int jday_first;
	jday_first= output->jday - 8;
#ifdef NRT
	jday_first = output->jday - 15;
#endif
	int nday;

	Obs *pobs = obs_list;
	while(pobs != NULL){
		if(!pobs->is_fill[band] && !pobs->is_outlier[band]){
			nday = pobs->jday - jday_first;		
			output->qa[col].valid_obs[band] = output->qa[col].valid_obs[band] | (1<<nday);

			if(output->jday == pobs->jday){
				output->qa[col].valid_obs_cur[band] = 1;
			}
		}

		pobs = pobs->next_v;
	}
}

void assess_band_qa(Output *output, int col, int band, int nobs)
{
	double rmse = output->rmse[col].band[band];
	double wd_r = output->wd_ref45[col].band[band];
	double wd_a = output->wd_wsa[col].band[band];

	int flag_band_qa;

	if(output->brdf[col].inv_status[band] == FULL_INVERSION){
		if(rmse < BAND_GOOD_RMSE){
			if(wd_r < BAND_NBAR_GOOD_WD){
				if(wd_a < BAND_WSA_GOOD_WD){
					flag_band_qa = 0;
				}
				else{
					flag_band_qa = 1;
				}
			}
			else{
				if(wd_a < BAND_WSA_GOOD_WD){
					flag_band_qa = 2;
				}
				else{
					flag_band_qa = 3;
				}
			}
		}
		else{
			if(wd_r < BAND_NBAR_GOOD_WD){
				if(wd_a < BAND_WSA_GOOD_WD){
					flag_band_qa = 4;
				}
				else{
					flag_band_qa = 5;
				}
			}
			else{
				if(wd_a < BAND_WSA_GOOD_WD){
					flag_band_qa = 6;
				}
				else{
					flag_band_qa = 7;
				}
			}
		}
	}
	else if(output->brdf[col].inv_status[band] == MAGN_INVERSION){
		if(nobs >= 7){
			flag_band_qa = 8;
		}
		else if(nobs > 0){
			flag_band_qa = 9;
		}
		else{
			flag_band_qa = QA_FILL;
		}
	}
	else{
		flag_band_qa = QA_FILL;
	}

	//
	if(flag_band_qa == 0 || flag_band_qa == 1 || flag_band_qa == 2 || flag_band_qa == 4){
		output->qa[col].band_qa[band] = 0;
	}
	else if(flag_band_qa == 3 || flag_band_qa == 5 || flag_band_qa == 6 || flag_band_qa == 7){
		output->qa[col].band_qa[band] = 1;
	}
	else if (flag_band_qa == 8){
		output->qa[col].band_qa[band] = 2;
	}
	else if (flag_band_qa != QA_FILL){
		output->qa[col].band_qa[band] = 3;
	}
	else{
		output->qa[col].band_qa[band] = QA_FILL;
	}

	// adjust by poor szn, quality terms, and valid_obs_curr
	if(output->qa[col].band_qa[band] == 0){
		if(output->qa[col].has_high_szn_obs   ||
		   output->qa[col].quality_term[band] ||
			 !output->qa[col].valid_obs_cur[band]){

			output->qa[col].band_qa[band] = 1;

		}
	}

	// mandatory qa
	if(output->brdf[col].inv_status[band] == FULL_INVERSION){
		output->qa[col].mandatory_qa[band] = 0;
	}
	else if(output->brdf[col].inv_status[band] == MAGN_INVERSION){
		output->qa[col].mandatory_qa[band] = 1;
	}
	else{
		output->qa[col].mandatory_qa[band] = QA_FILL;
	}

}

void assess_overall_qa(Output *output, int col, int is_snow)
{
	int band;
	
	// overall qa
	// overall qa not used
	/*double tot_rmse;
	double tot_wd_r;
	double tot_wd_a;
	int n_full;
	int n_magn;
	int n_fill;

	int flag_rmse = 100;
	int flag_wd_r = 100;
	int flag_wd_a = 100;

	int overall_qa;

	tot_rmse = 0.0;
	tot_wd_r = 0.0;
	tot_wd_a = 0.0;

	n_full = 0;
	n_magn = 0;
	n_fill = 0;

	for(band=0; band<NBAND; band++){
		if(output->brdf[col].inv_status[band] == FULL_INVERSION){
			tot_rmse += output->rmse[col].band[band]*output->rmse[col].band[band];
			tot_wd_r += output->wd_ref45[col].band[band]*output->wd_ref45[col].band[band];
			tot_wd_a += output->wd_bsa45[col].band[band]*output->wd_bsa45[col].band[band];
			n_full++;
		}
		else if(output->brdf[col].inv_status[band] == MAGN_INVERSION){
			n_magn++;
		}
		else{
			n_fill++;
		}
	}	

	if(n_full > 0){
		tot_rmse = sqrt(tot_rmse/(double)n_full);
		tot_wd_r = sqrt(tot_wd_r/(double)n_full);
		tot_wd_a = sqrt(tot_wd_a/(double)n_full);

		if(tot_rmse < MAX_GOOD_RMSE){
			flag_rmse = 0;
		}
		else if(tot_rmse < MAX_MODERATE_RMSE){
			flag_rmse = 1;
		}
		else if(tot_rmse < MAX_POOR_RMSE){
			flag_rmse = 2;
		}
		else{
			flag_rmse = 3;
		}

		if(tot_wd_r < 0 || tot_wd_r >= MAX_POOR_WD){
			flag_wd_r = 3;
		}
		else if(tot_wd_r < MAX_GOOD_WD){
			flag_wd_r = 0;
		}
		else if(tot_wd_r < MAX_MODERATE_WD){
			flag_wd_r = 1;
		}
		else if(tot_wd_r < MAX_POOR_WD){
			flag_wd_r = 2;
		}

		if(tot_wd_a < 0 || tot_wd_a >= MAX_POOR_WD){
			flag_wd_a = 3;
		}
		else if(tot_wd_a < MAX_GOOD_WD){
			flag_wd_a = 0;
		}
		else if(tot_wd_a < MAX_MODERATE_WD){
			flag_wd_a = 1;
		}
		else if(tot_wd_a < MAX_POOR_WD){
			flag_wd_a = 2;
		}

	}
	else{
		flag_rmse = 3;
		flag_wd_r = 3;
		flag_wd_a = 3;
	}

	overall_qa = (flag_rmse + flag_wd_r + flag_wd_a) / 3;

	if(n_full == NBAND && (overall_qa == 0 || overall_qa == 1)){
		overall_qa = 0;
	}	
	else if(n_fill == NBAND){
		overall_qa = QA_FILL;
	}
	else{
		overall_qa = 1;
	}

	output->qa[col].overall_qa = (int8)overall_qa;
	*/

	// broad band mandatory qa
	int i, j;
	int flags[3] = {0, 0, 0};

	for(i=0; i<3; i++){
		if(output->brdf[col].iso[NBAND+i] < 0.0){
			output->qa[col].mandatory_qa[NBAND+i] = QA_FILL;
			flags[i] = 1;
		}		
	}

	// sqs fixed a bug for snow mandatory qa for shortwave 3/6/2017
	int check_n[3] = {3, 4, 7};
	int check_i[3][7] = {{0,2,3,-1,-1,-1,-1},
											 {1,4,6,5,-1,-1,-1},
											 {0,1,2,4,6,3,5}};
	int flag;
	int flag_no_b6 = 0;
	if(is_snow){
		check_n[1] = 3;
		check_n[2] = 5;
	}
	else if(output->qa[col].mandatory_qa[5] == QA_FILL){
		check_n[1] = 3;
		check_n[2] = 6;
		flag_no_b6 = 1;
	}

	for(i=0; i<3; i++){
		if(flags[i]){	// already filled
			continue;
		}

		flag = 0;
		for(j=0; j<check_n[i]; j++){
			band = check_i[i][j];
			if(output->qa[col].band_qa[band] != 0 &&
				 output->qa[col].band_qa[band] != 1 ){
				flag = 1;
				break;	
			}	
		}
		if(i>0 && flag_no_b6){
			if(flag == 0){
				output->qa[col].mandatory_qa[NBAND+i] = 2;
			}
			else{
				output->qa[col].mandatory_qa[NBAND+i] = 3;
			}
		}
		else{
			if(flag == 0){
				output->qa[col].mandatory_qa[NBAND+i] = 0;
			}
			else{
				output->qa[col].mandatory_qa[NBAND+i] = 1;
			}
		}
	}

}

void n2b(float *pars, int8 *inv_status, int is_snow, float *vis, float *nir, float *sw, int has_intercept)
{
	static double coef[3][8] ={
	    { 0.3265, 0.0000, 0.4364, 0.2366, 0.0000, 0.0000, 0.0000, -0.0019 },
	    { 0.0000, 0.5447, 0.0000, 0.0000, 0.1363, 0.0469, 0.2536, -0.0068 },
	    { 0.3973, 0.2382, 0.3489, -0.2655, 0.1604, -0.0138, 0.0682, 0.0036 }};
	
	static double coef_no_b6[3][8] ={
	    { 0.3265, 0.0000, 0.4364, 0.2366, 0.0000, 0.0000, 0.0000, -0.0019 },
			{ 0.0000, 0.4597, 0.0000, 0.0000, 0.3398, 0.0000, 0.1968, -0.0092},
			{ 0.1806, 0.2062, 0.2327, 0.0834, 0.1770, 0.0000, 0.0928, -0.0054}};
	
	static double coef_snow[3][8] ={
			{ 0.3591, 0.0000, 0.5100, 0.1322, 0.0000, 0.0000, 0.0000, -0.0009},
			{ 0.0000, 0.6026, 0.0000, 0.0000, 0.2356, 0.0000, 0.1548, 0.0000},
			{ 0.1574, 0.2789, 0.3829, 0.0000, 0.1131, 0.0000, 0.0694, -0.0093}};	

	float *outp[3] = {vis, nir, sw};

	double wt_sum;
	double wt_tot;
	int i, b;
	double *pcoef;
	float ret;

	for(i=0; i<3; i++){
		wt_sum = 0.0;
		wt_tot = 0.0;
		ret = 0.0;

		if(is_snow){
			pcoef = coef_snow[i];
		}
		else{
			if(inv_status[5] == FILL_INVERSION){
				pcoef = coef_no_b6[i];
			}
			else{
				pcoef = coef[i];
			}
		}

		for(b=0; b<NBAND; b++){
			wt_tot += pcoef[b];
			if(inv_status[b] != FILL_INVERSION){
				wt_sum += pcoef[b];
				ret += pars[b] * pcoef[b];	
			}
		}	
		if(has_intercept){
			wt_tot += pcoef[b];
			wt_sum += pcoef[b];
			ret += pcoef[b];	
		}

		if(fabs(wt_sum - wt_tot) < 0.0001){
			valid_par(&ret, 0);
		}
		else{
			ret = DBL_FILL;
		}

		*(outp[i]) = ret;
	}	
}

void broad_brdf(Brdf *brdf, int is_snow)
{
	float *pars[3] = {brdf->iso, brdf->vol, brdf->geo};

	int i, j;
	int has_intercept;
	for(i=0; i<3; i++){
		if(i==0){
			has_intercept = 1;
		}
		else{
			has_intercept = 0;
		}

		n2b(pars[i], brdf->inv_status, is_snow, &(pars[i][NBAND]), &(pars[i][NBAND+1]), &(pars[i][NBAND+2]), has_intercept);

#ifdef DEBUG
	printf("ENTER BRD_BRDF PAR %d: VIS=%.4f, NIR=%.4f, SW=%.4f\n", i, pars[i][NBAND], pars[i][NBAND+1], pars[i][NBAND+2]);
#endif
	}

	// check
	int flag;
	for(i=0; i<3; i++){	// for vis, nir, sw
		flag = 0;
		for(j=0; j<3; j++){	// for iso, vol, geo
			if(pars[j][NBAND+i] < 0.0){
				flag = 1;
				break;
			}
		}
		if(flag){
			for(j=0; j<3; j++){	// for iso, vol, geo
				pars[j][NBAND+i] = DBL_FILL;	
			}	
			brdf->is_fill[NBAND+i] = 1;
			brdf->inv_status[NBAND+i] = FILL_INVERSION;
		}
		else{
			brdf->is_fill[NBAND+i] = 0;
		}
	}	

}

int valid_par(float *val, int limit_up)
{
	float v = *val;

	if(limit_up){
		if(v > 1.0){
			*val = DBL_FILL;
			return 0;
		}
	}

	if(v < -0.03){
		*val = DBL_FILL;
		return 0;
	}

	if(v < 0.0){
		v = 0.0;
	}

	*val = v;

	return 1;
}

void calc_albedo(Brdf *brdf, Albedo *albedo, Qa *qa, double lsn)
{
	int band;
	int flag = 0;

#ifdef NO_BROADBAND
	int nband = NBAND;
#else
	int nband = NBAND+3;
#endif

	for(band=0; band<nband; band++){
		//if(brdf->inv_status[band] == FILL_INVERSION){
		if(brdf->is_fill[band]){
			albedo->wsa[band] = DBL_FILL;			
			albedo->bsa[band] = DBL_FILL;			
			albedo->nbar[band] = DBL_FILL;			
		}
		else{
			float pars[3] = {brdf->iso[band], brdf->vol[band], brdf->geo[band]};

			albedo->wsa[band] = calc_wsa(pars);	
			albedo->bsa[band] = calc_bsa(pars, lsn);	
			albedo->nbar[band] = calc_nbar(pars, lsn);	

#ifdef DEBUG
				printf("ENTER CALC_ALB BAND %d BRDF=%.4f %.4f %.4f WSA=%.4f, BSA=%.4f, NBAR=%.4f, LSN=%.4f\n", 
								band, pars[0], pars[1], pars[2], albedo->wsa[band], albedo->bsa[band], albedo->nbar[band], lsn);
#endif

			if(!valid_par(&(albedo->wsa[band]), 1) ||
				 !valid_par(&(albedo->bsa[band]), 1) ||
				 !valid_par(&(albedo->nbar[band]), 1)	||
				 (int)lsn > 85 ){

				if(band >= NBAND){
					flag = 1;
					break;
				}

				brdf->is_fill[band] = 1;
				brdf->inv_status[band] = FILL_INVERSION;
				albedo->wsa[band] = DBL_FILL;
				albedo->bsa[band] = DBL_FILL;
				albedo->nbar[band] = DBL_FILL;
				qa->mandatory_qa[band] = QA_FILL;
				if(band<NBAND){
					qa->band_qa[band] = QA_FILL;
				}
			}
		
		}

	}	

	// if any broadband albedo is invalid, the quality of all the bands may not be good
	if(flag){
		for(band=0; band<nband; band++){
			brdf->is_fill[band] = 1;
			brdf->inv_status[band] = FILL_INVERSION;
			albedo->wsa[band] = DBL_FILL;
			albedo->bsa[band] = DBL_FILL;
			albedo->nbar[band] = DBL_FILL;
			qa->mandatory_qa[band] = QA_FILL;
			if(band<NBAND){
				qa->band_qa[band] = QA_FILL;
			}
		}
	}
}

void calc_lsn(Output *output, int jday)
{
	int year, doy;
	julday2doy(jday, &year, &doy);

	double lat_per_row = (North_boundary_lat - South_boundary_lat)/(double)output->ncol;	// should be nrow, but they're equal
	double lat = North_boundary_lat - lat_per_row * ((double)output->irow+0.5);

	output->lsn = get_local_szn_deg(lat, doy, 12.0);	
}

void platform_flag(Obs *obs_list, Output *output, int col)
{
	Obs *pobs = obs_list;

	int flag_am = 0;
	int flag_pm = 0;
	int flag = QA_FILL;

	while(pobs != NULL){
		if(pobs->platform == MODIS_TERRA){
			flag_am = 1;
		}
		if(pobs->platform == MODIS_AQUA){
			flag_pm = 1;
		}

		if(flag_am && flag_pm){
			break;
		}
			
		pobs = pobs->next_v;
	}

	if(flag_am && flag_pm){
		flag = TERRA_AQUA;	
	}
	else if(flag_am){
		flag = TERRA_ONLY;
	}
	else if(flag_pm){
		flag = AQUA_ONLY;
	}

	output->qa[col].platform_flag = flag;
}

void land_water_mask(Obs *obs_list, Output *output, int col)
{
	Obs *pobs = obs_list;

	int cnt[8] = {0,0,0,0,0,0,0,0};
	int mask = QA_FILL;

	while(pobs != NULL){
		mask = pobs->land_water_mask;
		if(mask < 8 && mask >= 0){
			cnt[mask] ++;
		}
			
		pobs = pobs->next_v;
	}

	int max = 0;

	int i;
	for(i=0; i<8; i++){
		if(cnt[i] > max){
			max = cnt[i];
			mask = i;
		}
	}

	if(max == 0){
		output->qa[col].land_water_mask = QA_FILL;
	}
	else{
		output->qa[col].land_water_mask = mask;
	}
}

void has_high_szn_obs(Obs *obs_list, Output *output, int col)
{
	Obs *pobs = obs_list;
	while(pobs != NULL){
		if(pobs->szn >= SZN_GOOD && pobs->szn < SZN_POOR){
			output->qa[col].has_high_szn_obs = 1;
			break;
		}
			
		pobs = pobs->next_v;
	}

}

void UpdateBKDB(Hdffile *backup_db, Output *output, Row *row)
{
	int col, band;	
	int season;

	for(col=0; col<row->ncol; col++){
		for(band=0; band<NBAND; band++){
			if(output->qa[col].band_qa[band] == 0 ||
		     output->qa[col].band_qa[band] == 1){

				if(!output->brdf[col].is_fill[band] &&
					 !EQ(output->brdf[col].iso[band], 0) &&
				   !EQ(output->brdf[col].vol[band], 0) &&
				   !EQ(output->brdf[col].geo[band], 0)){

					season = row->cells[col].season;
					assert(season >= 0 && season < 3);

					Brdf *pb[3] = {row->cells[col].brdf_bk->bk_s1, row->cells[col].brdf_bk->bk_s2,
									row->cells[col].brdf_bk->bk_s3};

					pb[season]->iso[band] = output->brdf[col].iso[band];
					pb[season]->vol[band] = output->brdf[col].vol[band];
					pb[season]->geo[band] = output->brdf[col].geo[band];
					pb[season]->is_fill[band] = 0;

				}
			}
		}
	}

	// moved to main()
	//assert(0 == WriteBKDBOneRow(backup_db, row));
}

void do_count(Row *input, Output *output)
{
	int col;
	int band;

	Counter *p = &(output->cnt);
	for(col=0; col<output->ncol; col++){
		if(!output->is_fill[col]){
			p->n_processed_pixel ++;
		}

		if(output->qa[col].land_water_mask == LAND){
			p->n_land_pixel ++;
		}
		else{
			p->n_notland_pixel ++;
		}

		if(input->cells[col].nobs	< 1){
			p->n_notprocessed_other += NBAND;
			continue;
		}

		p->n_total_obs += input->cells[col].nobs;

		for(band=0; band<NBAND; band++){
			if(output->qa[col].mandatory_qa[band] == 0){
				p->n_high_q ++;	
			}
			else if(output->qa[col].mandatory_qa[band] == 1){
				p->n_other_q ++;
			}
			else{
				p->n_notprocessed_cloud ++;
			}

			if(output->qa[col].band_qa[band] == 0 || 
				 output->qa[col].band_qa[band] == 1){
				p->n_new_brdf ++;
			}
			else if(output->qa[col].band_qa[band] == 2 || 
				 		  output->qa[col].band_qa[band] == 3){
				p->n_fixed_brdf ++;
			}
		}

	}	
}

int calc_ndsi(Obs *obs, double *ndsi)
{
	if(obs->is_fill[3] || obs->is_fill[5]){
		return -1;
	}

	*ndsi = (obs->ref[3] - obs->ref[5])/(obs->ref[3] + obs->ref[5]);

	return 0;
}

/* trim the observations by regions 
 * for polar tiles */
/* not used */
/*
#define MAX_OBS_PER_REGION 100
#define MAX_REGION 576
#define KEEP_OBS_PER_REGION 2

typedef struct{
	int azimuth;
	int zenith;
	int nobs;
	Obs *pobs[MAX_OBS_PER_REGION];
} Region;

Region regions[MAX_REGION];
int n_region = 0;

void init_regions()
{
	int i;
	n_region = 0;
	for(i=0; i<MAX_REGION; i++){
		regions[i].nobs = 0;
		regions[i].azimuth = -1;
		regions[i].zenith = -1;
	}
}

void add_to_regions(Obs *obs)
{
	if(!obs->selected){
		return;
	}

	if(obs->van < -180.0 || obs->van > 180.0){
		return;
	}
	if(obs->vzn < 0.0 || obs->vzn > 90.0){
		return;
	}

	int a = (int)(obs->van+180.0)/10;
	int z = (int)obs->vzn/5;

	assert(n_region < MAX_REGION);

	int i;
	for(i=0; i<n_region; i++){
		if(regions[i].azimuth == a && regions[i].zenith == z){
			int p = regions[i].nobs;
			assert(p < MAX_OBS_PER_REGION);

			regions[i].pobs[p] = obs;
			regions[i].nobs ++;
			return;
		}
	}

	// not found, new region
	regions[n_region].pobs[0] = obs;
	regions[n_region].nobs = 1;
	regions[n_region].azimuth = a;
	regions[n_region].zenith = z;
	//printf("adding region %d: %d %d\n", n_region, a, z);
	n_region ++;
}

void do_trim_obs()
{
	int i, j, k, b;
	int keep[KEEP_OBS_PER_REGION];
	double max_wt;

	for(i=0; i<n_region; i++){
		if(regions[i].nobs > KEEP_OBS_PER_REGION){
			for(j=0; j<KEEP_OBS_PER_REGION; j++){
				max_wt = 0.0;
				for(k=0; k<regions[i].nobs; k++){
					if(regions[i].pobs[k]->weight_v[0] > max_wt){
						//bug here!!! but if keep_obs_per_region == 2 it's okay
						if(j>0 && k == keep[j-1]){
							continue;
						}
						max_wt = regions[i].pobs[k]->weight_v[0];
						keep[j] = k;
					}
				}
			}
			
			for(k=0; k<regions[i].nobs; k++){
				int remove = 1;
				for(j=0; j<KEEP_OBS_PER_REGION; j++){
					if(k == keep[j]){
						remove = 0;	
						break;
					}
				}
				if(!remove){
					continue;
				}
				
				for(b=0; b<NBAND; b++){
					regions[i].pobs[k]->is_fill[b] = 1;
				}
			}

		}
	}
}

void do_trim_obs2()
{
	int i, b, m, n;
	Obs *ptmp;

	for(i=0; i<n_region; i++){
		if(regions[i].nobs > KEEP_OBS_PER_REGION){
			// sort by weight_v band 1
			for(m=0; m<regions[i].nobs; m++){
				for(n=m+1; n<regions[i].nobs; n++){
					if(regions[i].pobs[n]->weight_v[0] > regions[i].pobs[m]->weight_v[0]){
						ptmp = regions[i].pobs[m];
						regions[i].pobs[m] = regions[i].pobs[n];
						regions[i].pobs[n] = ptmp;
					}
				}
			}

				//remove
			for(b=0; b<NBAND; b++){
				for(m=KEEP_OBS_PER_REGION; m<regions[i].nobs; m++){
					regions[i].pobs[m]->is_fill[b] = 1;
					//printf("Removing region %d (azimuth: %d, zenith: %d) obs: ref band %d =%f, vaa=%f, vza=%f, wt=%f\n",
					//			 i, regions[i].azimuth*10, regions[i].zenith*5, b, regions[i].pobs[m]->ref[b], 
					//			 regions[i].pobs[m]->van, regions[i].pobs[m]->vzn, regions[i].pobs[m]->weight_v[b]);
				}
			}
		}
	}
}

void trim_obs(Obs *obs_list)
{
	Obs *pobs;

	init_regions();

	pobs = obs_list;
	while(pobs != NULL){
		add_to_regions(pobs);
		pobs = pobs->next_v;
	}

	do_trim_obs();
}
*/
