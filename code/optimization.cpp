#include "optimization.hpp"

// using namespace std;


double my_f(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    enew[i].length = exp(gsl_vector_get(v, i));
  }
  evo_tree new_tree(tree->nleaf, enew);
  if(lnl_type.model == MK){
      new_tree.mu = tree->mu;
  }else{
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}

double my_f_mu(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current branch lengths
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    enew[i].length = exp(gsl_vector_get(v, i));
  }
  evo_tree new_tree(tree->nleaf, enew);
  if(lnl_type.model == MK){
      new_tree.mu = exp(gsl_vector_get(v,tree->edges.size()-1));  // 0 to nedge-2 are epars, nedge - 1  is mu
  }else{
      new_tree.dup_rate = exp(gsl_vector_get(v,tree->edges.size()-1));
      new_tree.del_rate = exp(gsl_vector_get(v,tree->edges.size()));
      new_tree.chr_gain_rate = exp(gsl_vector_get(v,tree->edges.size() + 1));
      new_tree.chr_loss_rate = exp(gsl_vector_get(v,tree->edges.size()+2));
      new_tree.wgd_rate = exp(gsl_vector_get(v,tree->edges.size()+3));
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}


double my_f_cons(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }

  // parameters coming in are internal branch lengths followed by total time
  int count = 0;
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    if(enew[i].end > tree->nleaf - 1){
      enew[i].length = exp(gsl_vector_get(v, count));
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp(gsl_vector_get(v, count)));
  if(lnl_type.model == MK){
      new_tree.mu = tree->mu;
  }else{
      new_tree.dup_rate = tree->dup_rate;
      new_tree.del_rate = tree->del_rate;
      new_tree.chr_gain_rate = tree->chr_gain_rate;
      new_tree.chr_loss_rate = tree->chr_loss_rate;
      new_tree.wgd_rate = tree->wgd_rate;
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}

double my_f_cons_mu(const gsl_vector *v, void *params){
  GSL_PARAM* gsl_param = (GSL_PARAM*) params;
  evo_tree* tree = &(gsl_param->rtree);
  LNL_TYPE lnl_type = gsl_param->lnl_type;
  // evo_tree *tree = (evo_tree*) params;

  //create a new tree with the current parameters observing timing constraints
  vector<edge> enew;
  for(int i = 0; i < (tree->edges.size()); ++i){
    enew.push_back(tree->edges[i]);
  }

  // parameters coming in are internal branch lengths followed by mutation rate
  int count = 0;
  for(int i = 0; i < (tree->edges.size())-1; ++i){
    if(enew[i].end > tree->nleaf - 1){
      enew[i].length = exp(gsl_vector_get(v, count));
      count++;
    }else{
      enew[i].length = 0;
    }
  }

  evo_tree new_tree(tree->nleaf, enew, exp(gsl_vector_get(v, count)));
  if(lnl_type.model == MK){
      new_tree.mu = exp(gsl_vector_get(v, count + 1));
  }
  else{
      new_tree.dup_rate = exp(gsl_vector_get(v, count + 1));
      new_tree.del_rate = exp(gsl_vector_get(v, count + 2));
      new_tree.chr_gain_rate = exp(gsl_vector_get(v, count + 3));
      new_tree.chr_loss_rate = exp(gsl_vector_get(v, count + 4));
      new_tree.wgd_rate = exp(gsl_vector_get(v, count + 5));
  }

  return -1.0 * get_likelihood_revised(new_tree, gsl_param->vobs, lnl_type);
}


// given a tree, maximise the branch lengths(and optionally mu) assuming branch lengths are independent or constrained in time
// use GSL simplex optimization
void max_likelihood(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, const vector<double>& tobs, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double& min_nlnl, const double& ssize){
  int debug = 0;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  /* save original handler, install new handler */
  gsl_error_handler_t  *old_handler = gsl_set_error_handler(&my_err_handler);

  int model = lnl_type.model;
  int cn_max = lnl_type.cn_max;
  int only_seg = lnl_type.only_seg;
  int is_total = lnl_type.is_total;

  int cons = lnl_type.cons;
  int maxj = opt_type.maxj;
  double tolerance = opt_type.tolerance;
  int miter = opt_type.miter;

  int npar_ne;
  int npar;
  gsl_vector *x;

  int nedge = 2 * rtree.nleaf - 2;
  int nintedge = rtree.nleaf - 2;

  if(!cons){
    npar_ne = nedge - 1 ;
    if(!maxj){
      npar = npar_ne;
    }else{
      if(model == MK){
          npar = npar_ne + 1;
      }
      if(model == BOUNDT){
          npar = npar_ne + 5;
      }
    }

    // initialise the best guess branch length and mu if required
    x = gsl_vector_alloc(npar);
    for(int i = 0; i < npar_ne; ++i){
      gsl_vector_set(x, i, log(rtree.edges[i].length));
    }
    if(maxj){
      if(model == MK){
          gsl_vector_set(x, npar_ne, log(rtree.mu));
      }
      if(model == BOUNDT){
          gsl_vector_set(x, npar_ne, log(rtree.dup_rate));
          gsl_vector_set(x, npar_ne + 1, log(rtree.del_rate));
          gsl_vector_set(x, npar_ne + 2, log(rtree.chr_gain_rate));
          gsl_vector_set(x, npar_ne + 3, log(rtree.chr_loss_rate));
          gsl_vector_set(x, npar_ne + 4, log(rtree.wgd_rate));
      }
      minex_func.f = my_f_mu;
    }else{
      minex_func.f = my_f;
    }
  }else{
    npar_ne = nintedge + 1;
    if(!maxj){
      npar = npar_ne;
    }else{
        if(model == MK){
            npar = npar_ne + 1;
        }
        if(model == BOUNDT){
            npar = npar_ne + 5;
        }
    }

    x = gsl_vector_alloc(npar);
    // initialise with internal edges
    vector<edge*> intedges = rtree.get_internal_edges();
    for(int i = 0; i < nintedge; ++i){
      gsl_vector_set(x, i, log(intedges[i]->length));
    }

    // initialise with current total tree time
    gsl_vector_set(x, nintedge, log(get_total_time(rtree.get_node_times(), lnl_type.max_tobs)));
    if(maxj){
      if(model == MK){
          gsl_vector_set(x, npar_ne, log(rtree.mu));
      }
      if(model == BOUNDT){
          gsl_vector_set(x, npar_ne, log(rtree.dup_rate));
          gsl_vector_set(x, npar_ne + 1, log(rtree.del_rate));
          gsl_vector_set(x, npar_ne + 2, log(rtree.chr_gain_rate));
          gsl_vector_set(x, npar_ne + 3, log(rtree.chr_loss_rate));
          gsl_vector_set(x, npar_ne + 4, log(rtree.wgd_rate));
      }
      minex_func.f = my_f_cons_mu;
    }else{
      minex_func.f = my_f_cons;
    }
  }

  // Set initial step sizes to 1
  //cout << "numbers: nedge / nintedge / npar / x->size: " << nedge << "\t" << nintedge << "\t" << npar << "\t" << x->size << endl;
  gsl_vector* ss = gsl_vector_alloc(npar);
  gsl_vector_set_all(ss, ssize);

  // Initialize method and iterate
  // evo_tree *p = &rtree;
  GSL_PARAM gsl_param = {rtree, vobs, lnl_type};
  GSL_PARAM* p = &gsl_param;
  void* pv = p;
  minex_func.n = npar;
  minex_func.params = pv;

  s = gsl_multimin_fminimizer_alloc(T, npar);
  gsl_multimin_fminimizer_set(s, &minex_func, x, ss);

  size_t iter = 0;
  int status;
  double size;
  do{
    iter++;
    status = gsl_multimin_fminimizer_iterate(s);
    if(status) break;

    size = gsl_multimin_fminimizer_size(s);
    status = gsl_multimin_test_size(size, tolerance);
    // gsl_set_error_handler_off();

    if(debug){
      printf("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
	      iter,
	      exp(gsl_vector_get(s->x, 0)),
	      exp(gsl_vector_get(s->x, s->x->size-1)),
	      s->fval, size);
    }

    if(status == GSL_SUCCESS){
      if(debug){
        	printf("converged to minimum at\n");

        	printf("%5lu, %10.3e, %10.3e, f() = %7.3f, size = %.3f\n",
        		iter,
        		exp(gsl_vector_get(s->x, 0)),
        		exp(gsl_vector_get(s->x, s->x->size-1)),
        		s->fval, size);
        }
    }
  }
  while(status == GSL_CONTINUE && iter < miter);

  if(status == GSL_CONTINUE ){
    cout << "### WARNING: maximum likelihood did not converge" << endl;
  }

  if(!cons){
    for(int i = 0; i < npar_ne; ++i){
      rtree.edges[i].length = exp(gsl_vector_get(s->x, i));
    }
    if(maxj){
        if(model == MK){
            rtree.mu = exp(gsl_vector_get(s->x, npar_ne));
        }
        if(model == BOUNDT){
            rtree.dup_rate = exp(gsl_vector_get(s->x, npar_ne));
            rtree.del_rate = exp(gsl_vector_get(s->x, npar_ne + 1));
            rtree.chr_gain_rate = exp(gsl_vector_get(s->x, npar_ne+2));
            rtree.chr_loss_rate = exp(gsl_vector_get(s->x, npar_ne+3));
            rtree.wgd_rate = exp(gsl_vector_get(s->x, npar_ne+4));
            rtree.mu = 0;
        }
    }

    min_nlnl = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    /* restore original handler */
    gsl_set_error_handler(old_handler);

  }else{
    int count = 0;
    // update internal lengths
    for(int i = 0; i < nedge - 1 ; ++i){
      if(rtree.edges[i].end > rtree.nleaf - 1){
        	rtree.edges[i].length = exp(gsl_vector_get(s->x, count));
        	count++;
      }else{
	        rtree.edges[i].length = 0;
      }
    }

    // Fill external edge lengths by looping over nodes
    double total_time = exp(gsl_vector_get(s->x, count));
    for(int i = 0; i < rtree.nleaf - 1; ++i){
      vector<int> es = rtree.get_ancestral_edges(rtree.nodes[i].id);
      reverse(es.begin(), es.end());

      rtree.edges[es.back()].length = total_time + tobs[ rtree.nodes[i].id ];
      for(int j = 0; j < es.size()-1; ++j){
        rtree.edges[es.back()].length -= rtree.edges[es[j]].length;
      }
    }

    if(maxj){
        if(model == MK){
            rtree.mu = exp(gsl_vector_get(s->x, npar_ne));
        }
        if(model == BOUNDT){
            rtree.dup_rate = exp(gsl_vector_get(s->x, npar_ne));
            rtree.del_rate = exp(gsl_vector_get(s->x, npar_ne + 1));
            rtree.chr_gain_rate = exp(gsl_vector_get(s->x, npar_ne+2));
            rtree.chr_loss_rate = exp(gsl_vector_get(s->x, npar_ne+3));
            rtree.wgd_rate = exp(gsl_vector_get(s->x, npar_ne+4));
            rtree.mu = 0;
        }
    }

    min_nlnl = s->fval;
    gsl_vector_free(x);
    gsl_vector_free(ss);
    gsl_multimin_fminimizer_free(s);
    /* restore original handler */
    gsl_set_error_handler(old_handler);

  }

}


/*****************************************************
    One dimensional optimization with Brent method
*****************************************************/

double computeFunction(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double value, int type){
    // Update the branch
    switch(type){
        case -1:{
            rtree.current_it->length = value;
            rtree.current_it_back->length = value;
            // Need to update node times and ages?
            rtree.edges[rtree.current_eid].length = value;
            break;
        }
        case 0:{
            rtree.mu = value;
            break;
        }
        case 1:{
            rtree.dup_rate = value;
            break;
        }
        case 2:{
            rtree.del_rate = value;
            break;
        }
        case 3:{
            rtree.chr_gain_rate = value;
            break;
        }
        case 4:{
            rtree.chr_loss_rate = value;
            break;
        }
        case 5:{
            rtree.wgd_rate = value;
            break;
        }
        default:
            break;
    }

    double nlnl = 0.0;
    if(lnl_type.model == DECOMP){
        nlnl = -get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
    }else{
        nlnl = -get_likelihood_revised(rtree, vobs, lnl_type);
    }
    // cout << "\n result of computeFunction " << nlnl << endl;

    return nlnl;
}


#define ITMAX 100
#define CGOLD 0.3819660
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d)(a)= (b);(b)= (c);(c)= (d);
#define SIGN(a,b)((b) >= 0.0 ? fabs(a) : -fabs(a))

/* Brents method in one dimension */
double brent_opt(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double ax, double bx, double cx, double tol,
                          double *foptx, double *f2optx, double fax, double fbx, double fcx){
	int iter;
	double a,b,d = 0,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double xw,wv,vx;
	double e = 0.0;

	a= (ax < cx ? ax : cx);
	b= (ax > cx ? ax : cx);
	x=bx;
	fx=fbx;
	if(fax < fcx){
		w=ax;
		fw=fax;
		v=cx;
		fv=fcx;
	} else{
		w=cx;
		fw=fcx;
		v=ax;
		fv=fax;
	}

	for(iter=1;iter<=ITMAX;iter++){
		xm = 0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if(fabs(x-xm) <= (tol2 - 0.5*(b-a))){
			*foptx = fx;
			xw = x-w;
			wv = w-v;
			vx = v-x;
			*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/
			         (v*v*xw + x*x*wv + w*w*vx);
			return x;
		}

		if(fabs(e) > tol1){
			r= (x-w)*(fx-fv);
			q= (x-v)*(fx-fw);
			p= (x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if(q > 0.0)
				p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if(fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e= (x >= xm ? a-x : b-x));
			else{
				d=p/q;
				u=x+d;
				if(u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else{
			d=CGOLD*(e= (x >= xm ? a-x : b-x));
		}

		u = (fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, u, type);
		if(fu <= fx){
			if(u >= x)
				a=x;
			else
				b=x;

			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else{
			if(u < x)
				a=u;
			else
				b=u;
			if(fu <= fw || w == x){
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else
				if(fu <= fv || v == x || v == w){
					v=u;
					fv=fu;
				}
		}
	}

	*foptx = fx;
	xw = x-w;
	wv = w-v;
	vx = v-x;
	*f2optx = 2.0*(fv*xw + fx*wv + fw*vx)/(v*v*xw + x*x*wv + w*w*vx);

	return x;
}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef SIGN
#undef GOLD
#undef GLIMIT
#undef TINY


double minimizeOneDimen(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double xmin, double xguess, double xmax, double tolerance, double *fx, double *f2x){
	double eps, optx, ax, bx, cx, fa, fb, fc;
	//int    converged;	/* not converged error flag */

	/* first attempt to bracketize minimum */
	if(xguess < xmin) xguess = xmin;
	if(xguess > xmax) xguess = xmax;
	eps = xguess*tolerance*50.0;
	ax = xguess - eps;
	if(ax < xmin) ax = xmin;
	bx = xguess;
	cx = xguess + eps;
	if(cx > xmax) cx = xmax;

	/* check if this works */
    // compute fb first to save some computation, if any
	fb = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, bx, type);
	fa = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, ax, type);
	fc = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, cx, type);

	/* if it works use these borders else be conservative */
	if((fa < fb) ||(fc < fb)){
		if(ax != xmin) fa = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, xmin, type);
		if(cx != xmax) fc = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, xmax, type);
		ax = xmin;
		cx = xmax;
	}
	/*
	const int MAX_ROUND = 10;
	for(i = 0;((fa < fb-tolerance) ||(fc < fb-tolerance)) &&(i < MAX_ROUND); i++){
		// brent method require that fb is smaller than both fa and fc
		// find some random values until fb achieve this
			bx = (((double)rand()) / RAND_MAX)*(cx-ax) + ax;
			fb = computeFunction(bx);
	}*/

/*
	if((fa < fb) ||(fc < fb)){
		if(fa < fc){ bx = ax; fb = fa; } else{ bx = cx; fb = fc; }
		//cout << "WARNING: Initial value for Brent method is set at bound " << bx << endl;
	}*/
	//	optx = brent_opt(xmin, xguess, xmax, tolerance, fx, f2x, fa, fb, fc);
	//} else
	optx = brent_opt(rtree, vobs, obs_decomp, comps, lnl_type, type, ax, bx, cx, tolerance, fx, f2x, fa, fb, fc);
  if(*fx > fb) // if worse, return initial value
  {
      *fx = computeFunction(rtree, vobs, obs_decomp, comps, lnl_type, bx, type);
      return bx;
  }

	return optx; /* return optimal x */
}


// return log likelihood
double optimize_mutation_rates(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double tolerance){
    int debug = 0;

    if(debug){
        cout << "\tUsing Brent method to optimize the likelihood of mutation rate" << endl;
    }
    double negative_lh = MAX_NLNL;
    double ferror, optx;

    if(lnl_type.model == MK){
        optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, 0, RATE_MIN, rtree.mu, RATE_MAX, tolerance, &negative_lh, &ferror);
    }else{
        optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, 1, RATE_MIN, rtree.dup_rate, RATE_MAX, tolerance, &negative_lh, &ferror);
        rtree.dup_rate = optx;
        if(debug){
            cout << "\tmax Brent logl: " << -negative_lh << " optimized duplication rate " << optx << endl;
        }

        optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, 2, RATE_MIN, rtree.del_rate, RATE_MAX, tolerance, &negative_lh, &ferror);
        rtree.del_rate = optx;
        if(debug){
            cout << "\tmax Brent logl: " << -negative_lh << " optimized deletion rate " << optx << endl;
        }

        if(!lnl_type.only_seg){
            optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, 3, RATE_MIN, rtree.chr_gain_rate, RATE_MAX, tolerance, &negative_lh, &ferror);
            rtree.chr_gain_rate = optx;
            if(debug){
                cout << "\tmax Brent logl: " << -negative_lh << " optimized chromosome gain rate " << optx << endl;
            }

            optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, 4, RATE_MIN, rtree.chr_loss_rate, RATE_MAX, tolerance, &negative_lh, &ferror);
            rtree.chr_loss_rate = optx;
            if(debug){
                cout << "\tmax Brent logl: " << -negative_lh << " optimized chromosome loss rate " << optx << endl;
            }

            optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, 5, RATE_MIN, rtree.wgd_rate, RATE_MAX, tolerance, &negative_lh, &ferror);
            rtree.wgd_rate = optx;
            if(debug){
                cout << "\tmax Brent logl: " << -negative_lh << " optimized WGD rate " << optx << endl;
            }
        }
    }
    return -negative_lh;
}


// return log likelihood
double optimize_one_branch(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, 
LNL_TYPE& lnl_type, double tolerance, Node* node1, Node* node2){
    int debug = 0;
    if(debug){
        cout << "\tOptimizing the branch " << node1->id + 1 << ", " << node2->id + 1 << endl;
    }

    assert(!((node1->id == rtree.nleaf && node2->id == rtree.nleaf - 1) || (node2->id == rtree.nleaf && node1->id == rtree.nleaf - 1)));
        
    rtree.current_it = (Neighbor*) node1->findNeighbor(node2);
    assert(rtree.current_it);
    rtree.current_it_back = (Neighbor*) node2->findNeighbor(node1);
    assert(rtree.current_it_back);
    assert(rtree.current_it->length == rtree.current_it_back->length);

    int eid = rtree.get_edge_id(node1->id, node2->id);
    assert(eid >= 0);
    rtree.current_eid = eid;

    double current_len = rtree.edges[eid].length;
    assert(current_len >= 0.0);
    double negative_lh = MAX_NLNL;
    double ferror, optx;

    // Brent method
    optx = minimizeOneDimen(rtree, vobs, obs_decomp, comps, lnl_type, -1, BLEN_MIN, current_len, BLEN_MAX, tolerance, &negative_lh, &ferror);

    rtree.current_it->length = optx;
    rtree.current_it_back->length = optx;
    rtree.edges[eid].length = optx;

    if(debug){
        cout << "\tUsing Brent method to optimize the likelihood of one branch length of edge " << eid + 1 << endl;
        cout << "\tmax Brent logl: " << -negative_lh << " optimized branch length " << optx << endl;
    }

    return -negative_lh;
}


void compute_best_traversal(evo_tree& rtree, NodeVector &nodes, NodeVector &nodes2){
    Node* farleaf = rtree.find_farthest_leaf();

    // double call to farthest leaf to find the longest path on the tree
    rtree.find_farthest_leaf(farleaf);

    rtree.get_preorder_branches(nodes, nodes2, farleaf);
}


double optimize_all_branches(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, 
LNL_TYPE& lnl_type, int my_iterations, double tolerance){
    int debug = 0;

    NodeVector nodes, nodes2;
    compute_best_traversal(rtree, nodes, nodes2);

    if(debug){
        cout << "\nOptimizing branch lengths (max " << my_iterations << " loops)..." << endl;
        cout << "nodes in best traversal: ";
        for(auto n : nodes){
            cout << "\t" << n->id + 1;
        }
        cout << endl;
        cout << "nodes2 in best traversal: ";
        for(auto n : nodes2){
            cout << "\t" << n->id + 1;
        }
        cout << endl;
    }

    int model = lnl_type.model;
    double tree_lh = rtree.score;

    // if(model == DECOMP){
    //     tree_lh = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
    // }else{
    //     tree_lh = get_likelihood_revised(rtree, vobs, lnl_type);
    // }
    if(debug){
        cout << "Initial tree log-likelihood: " << tree_lh << endl;
    }

    DoubleVector lenvec;
    double new_tree_lh = 0.0;
    for(int i = 0; i < my_iterations; i++){
    	  save_branch_lengths(rtree, lenvec, 0);

        for(int j = 0; j < nodes.size(); j++){
            if(debug){
                cout << " optimizing branch " << nodes[j]->id << " " << nodes2[j]->id << endl;
            }   
            // skip normal branch, which will always be 0
            if(nodes[j]->id == rtree.nleaf - 1 || nodes2[j]->id == rtree.nleaf - 1)   continue;    
            new_tree_lh = optimize_one_branch(rtree, vobs, obs_decomp, comps, lnl_type, tolerance, (Node* )nodes[j],(Node* )nodes2[j]);
        }
        
        if(debug){
            cout << "tree log likelihood after iteration " << i + 1 << " : " << new_tree_lh << endl;
        }

        if(new_tree_lh < tree_lh - tolerance * 0.1){
            // IN RARE CASE: tree log-likelihood decreases, revert the branch length and stop
            if(debug){
                cout << "tree log-likelihood decreases" << endl;
                cout << "NOTE: Restoring branch lengths as tree log-likelihood decreases after branch length optimization: " << tree_lh << " -> " << new_tree_lh << endl;
            }
            
            restore_branch_lengths(rtree, lenvec);

            double max_delta_lh = 1.0;

            if(model == DECOMP){
              new_tree_lh = get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
            }else{
              new_tree_lh = get_likelihood_revised(rtree, vobs, lnl_type);
            }

            if(fabs(new_tree_lh - tree_lh) > max_delta_lh){
              cout << endl;
              cout << "new_tree_lh: " << new_tree_lh << "   tree_lh: " << tree_lh << endl;
            }
        	assert(fabs(new_tree_lh - tree_lh) < max_delta_lh);

        	return new_tree_lh;
        }

        // only return if the new_tree_lh >= tree_lh!(in rare case that likelihood decreases, continue the loop)
        if(tree_lh <= new_tree_lh && new_tree_lh <= tree_lh + tolerance){
          if(debug) cout << "tree log-likelihood increases" << endl;
        	return new_tree_lh;
        }

        tree_lh = new_tree_lh;
    }

    // recompute node times and ages to be consistent with branch lengths
    rtree.calculate_node_times();
    rtree.calculate_age_from_time();

    if(debug) cout << "current score " << tree_lh << endl;

    return tree_lh;
}


/*****************************************************
    L-BFGS-B method
*****************************************************/
void update_variables(evo_tree& rtree, int model, int cons, int maxj, double *x){
    int debug = 0;

    int nedge = 2 * rtree.nleaf - 2;
    // create a new tree from current value of parameters
    vector<edge> enew;
    for(int i = 0; i < nedge; ++i){
      enew.push_back(rtree.edges[i]);
    }

    if(!cons){
      // The first element of x is not used for optimization
      // The index of x is added by 1 compared with index for simplex method
      for(int i = 0; i < nedge - 1 ; ++i){
          enew[i].length = x[i + 1];
      }
      evo_tree new_tree(rtree.nleaf, enew);
      if(maxj){
          if(model == MK){
              new_tree.mu = x[nedge];
              if(debug){
                  for(int i = 0; i < nedge + 1; i++){ cout << x[i] << '\n';}
                  cout << "mu value so far: " << new_tree.mu << endl;
              }
          }else{
              new_tree.dup_rate = x[nedge];
              new_tree.del_rate = x[nedge + 1];
              new_tree.chr_gain_rate = x[nedge+2];
              new_tree.chr_loss_rate = x[nedge+3];
              new_tree.wgd_rate = x[nedge+4];
              new_tree.mu = 0;
              if(debug){
                  for(int i = 0; i < nedge+2; i++){ cout << x[i] << '\n';}
                  cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                  cout << "del_rate value so far: " << new_tree.del_rate << endl;
                  cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                  cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                  cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
              }
          }
      }
      rtree = evo_tree(new_tree);
    }else{
      int count = 0;
      for(int i = 0; i < nedge - 1 ; ++i){
        if(enew[i].end > rtree.nleaf - 1){
          	enew[i].length = x[count + 1];
          	count++;
        }else{
  	        enew[i].length = 0;
        }
      }
      if(debug){
          cout << "total height so far: " << x[count + 1] << endl;
      }
      evo_tree new_tree(rtree.nleaf, enew, x[count + 1]);
      if(maxj){
        int nintedge = rtree.nleaf - 2;

        if(model == MK){
            new_tree.mu = x[nintedge + 2];
            if(debug){
                for(int i = 0; i <= nintedge + 2; i++){ cout << x[i] << '\n';}
                cout << "mu value so far: " << new_tree.mu << endl;
            }
        }
        else{
            new_tree.dup_rate = x[nintedge + 2];
            new_tree.del_rate = x[nintedge + 3];
            new_tree.chr_gain_rate = x[nintedge + 4];
            new_tree.chr_loss_rate = x[nintedge + 5];
            new_tree.wgd_rate = x[nintedge + 6];
            new_tree.mu = 0;
            if(debug){
                for(int i = 0; i <= nintedge + 6; i++){ cout << x[i] << '\n';}
                cout << "dup_rate value so far: " << new_tree.dup_rate << endl;
                cout << "del_rate value so far: " << new_tree.del_rate << endl;
                cout << "chr_gain_rate value so far: " << new_tree.chr_gain_rate << endl;
                cout << "chr_loss_rate value so far: " << new_tree.chr_loss_rate << endl;
                cout << "wgd_rate value so far: " << new_tree.wgd_rate << endl;
            }
         }
      }
      rtree = evo_tree(new_tree);
    }
}


// Update the tree after each iteration in the BFGS optimization
// Estimate ratios rather than branch length in order to avoid negative terminal branch lengths
// Sort node times in increasing order and take the first Ns intervals
void update_variables_transformed(evo_tree& rtree, double *x, LNL_TYPE& lnl_type, OPT_TYPE& opt_type){
    int debug = 0;
    if(debug){
        cout << "Update the tree after each iteration in the BFGS optimization" << endl;
        cout << "tree before optimization " << rtree.make_newick() << endl;
        cout << "all the node ages so far: ";
        for(auto n : rtree.nodes){
          cout << "\t" << n.age;
        }
        cout << endl;
    }

    int cn_max = lnl_type.cn_max;
    int only_seg = lnl_type.only_seg;
    int is_total = lnl_type.is_total;

    int npar_ne = 0;
    if(opt_type.opt_one_branch){
        if(debug){
          cout << "transform only one branch " << rtree.current_eid + 1 << endl;
        }
        if(!lnl_type.cons){
            npar_ne = 1;
            rtree.edges[rtree.current_eid].length = x[1];
            rtree.current_it->length = x[1];
            rtree.current_it_back->length = x[1];
        }else{
            npar_ne = 1;
            vector<double> ratios = rtree.get_ratio_from_age();

            // ratios[0] = x[1]; // need to update root edge
            int eend = rtree.edges[rtree.current_eid].end;
            assert(eend > rtree.nleaf);
            int nid = eend - rtree.root_node_id;
            ratios[nid] = x[1];
            // update based on all ratios together to avoid inconsistencies
            rtree.update_edges_from_ratios(ratios, lnl_type.knodes);
            // rtree.update_edge_from_ratio(x[1], rtree.current_eid);
            // get length of optimized edge
            double blen = rtree.edges[rtree.current_eid].length;
            rtree.current_it->length = blen;
            rtree.current_it_back->length = blen;

            if(debug){
              cout << "estimated x: ";
              for(int i = 0; i < npar_ne; i++){
                  double val = x[i + 1];
                  cout << "\t" << val;
              }
              cout << endl;
              cout << "tree after optimization " << rtree.make_newick() << endl;
              cout << "all the node ages so far: ";
              for(auto n : rtree.nodes){
                cout << "\t" << n.age;
              }
              cout << endl;
            }
        }
    }else{
        if(!lnl_type.cons){
            // npar_ne = nedge - 1;
            npar_ne = 2 * rtree.nleaf - 3;
            // The first element of x is not used for optimization
            // The index of x is added by 1 compared with index for simplex method
            for(int i = 0; i < npar_ne; i++){
              rtree.edges[i].length = x[i + 1];
            }
        }else{
            // npar_ne = nintedge + 1;
            npar_ne = rtree.nleaf - 1;

            // store original ratios
            vector<double> ratios = rtree.get_ratio_from_age();
           
            double min_root = lnl_type.max_tobs + rtree.nleaf * BLEN_MIN;

            if(debug){
                int nratios = rtree.nleaf - 1;
                cout << "There are " << nratios << " ratios" << endl;
                cout << "Original values of estimated variables: " << endl;
                for(int i = 0; i < nratios; i++){
                    cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
                }
                cout << "estimated x: ";
                for(int i = 0; i < npar_ne; i++){
                    double val = x[i + 1];
                    cout << "\t" << val;
                }
                cout << endl;
                cout << "min age of root allowed " << min_root << endl;
            }
            // The estimated value may be nan
            for(int i = 0; i < npar_ne; i++){
                double val = x[i + 1];
                if(std::isnan(val)){  // return previous values
                    if(debug) cout << "nan returned in BFGS!" << endl;
                    val = ratios[i];
                }
                ratios[i] = val;
            }

            // If the optimized root is close to the boundary, revert to original values to avoid very small branches
            if(fabs(x[1] - min_root) > SMALL_VAL){
              // cout << "\n\nupdate all branches" << endl;
              rtree.update_edges_from_ratios(ratios, lnl_type.knodes);
            }

            if(!is_age_time_consistent(rtree.get_node_times(), rtree.get_node_ages())){
              cout << "Updated tree has inconsistent time/ages after updating the optimized tree!" << endl;
              rtree.print();
              string newick = rtree.make_newick(8);
              cout << newick << endl;
            }

            if(debug){
                cout << "Current values of estimated variables: " << endl;
                for(int i = 0; i < ratios.size(); i++){
                    cout << i + 1 << "\t" << "\t" << ratios[i] << endl;
                }
            }
        }

        if(debug){
          cout << "New branch lengths: \n";
          for(int i = 0; i < rtree.edges.size(); i++){
              cout << i + 1 << "\t" << "\t" << rtree.edges[i].length << endl;
          }
        }
    }

    if(lnl_type.cons && !is_tip_age_valid(rtree.get_node_ages(), opt_type.tobs)){
      rtree.print();
      cout << rtree.make_newick() << endl;
      cout << "Tip timings inconsistent with observed data after updating the optimized tree!" << endl;
      exit(1);
    }

    if(opt_type.maxj){
        if(lnl_type.model == MK){
            // nintedge + 2 for constrained branches
            rtree.mu = x[npar_ne + 1];
            if(debug){
                for(int i = 0; i <= npar_ne + 1; i++){ cout << x[i] << '\n'; }
                cout << "mu value so far: " << rtree.mu << endl;
            }
        }else{
            rtree.dup_rate = x[npar_ne + 1];
            rtree.del_rate = x[npar_ne + 2];
            if(!only_seg){
                rtree.chr_gain_rate = x[npar_ne + 3];
                rtree.chr_loss_rate = x[npar_ne + 4];
                rtree.wgd_rate = x[npar_ne + 5];
            }

            if(debug){
                for(int i = 0; i <= npar_ne + 1; i++){ cout << x[i] << '\n';}
                cout << "dup_rate value so far: " << rtree.dup_rate << endl;
                cout << "del_rate value so far: " << rtree.del_rate << endl;
                if(!only_seg){
                    cout << "chr_gain_rate value so far: " << rtree.chr_gain_rate << endl;
                    cout << "chr_loss_rate value so far: " << rtree.chr_loss_rate << endl;
                    cout << "wgd_rate value so far: " << rtree.wgd_rate << endl;
                }
            }
        }
    }
}



/**
    the target function which needs to be optimized (negative log likelihood)
    @param x the input vector x
    @return the function value at x (negative log likelihood)
*/
double targetFunk(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, 
LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double x[]){
    update_variables_transformed(rtree, x, lnl_type, opt_type);

    if(lnl_type.model == DECOMP){
        return -1.0 * get_likelihood_decomp(rtree, vobs, obs_decomp, comps, lnl_type);
    }else{
        return -1.0 * get_likelihood_revised(rtree, vobs, lnl_type);
    }
}


/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
double derivativeFunk(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, 
LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int ndim, double x[], double dfx[]){
  int debug = 0;

	double *h = new double[ndim + 1];
  double temp;
  int dim;

  double fx = targetFunk(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, x);

	for(dim = 1; dim <= ndim; dim++){
		temp = x[dim];
		h[dim] = ERROR_X * fabs(temp);
		if(h[dim] == 0.0) h[dim] = ERROR_X;
		x[dim] = temp + h[dim];
		h[dim] = x[dim] - temp;

    dfx[dim] = (targetFunk(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, x));
		x[dim] = temp;
	}

	for(dim = 1; dim <= ndim; dim++){
        dfx[dim] = (dfx[dim] - fx) / h[dim];
        if(debug){
            cout << "dfx[dim] " << dfx[dim] << endl;
        }
    }

  delete [] h;

	return fx;
}


double optimFunc(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *vars){
    return targetFunk(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, vars-1);
}

double optimGradient(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *x, double *dfx){
    return derivativeFunk(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, nvar, x - 1, dfx - 1);
}


void lbfgsb(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT){
	char task[60];
	double f, *g, dsave[29], *wa;
	int tr = -1, iter = 0, *iwa, isave[44], lsave[4];

	/* shut up gcc -Wall in 4.6.x */

	for(int i = 0; i < 4; i++) lsave[i] = 0;

	if(n == 0){ /* not handled in setulb */
		*fncount = 1;
		*grcount = 0;

    *Fmin = optimFunc(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, n, u);
		strcpy(msg, "NOTHING TO DO");
		*fail = 0;
		return;
	}

	if(nREPORT <= 0){
		cerr << "REPORT must be > 0(method = \"L-BFGS-B\")" << endl;
		exit(1);
	}

	switch(trace){
  	case 2: tr = 0; break;
  	case 3: tr = nREPORT; break;
  	case 4: tr = 99; break;
  	case 5: tr = 100; break;
  	case 6: tr = 101; break;
  	default: tr = -1; break;
	}

	*fail = 0;
	g = (double*) malloc(n * sizeof(double));
	/* this needs to be zeroed for snd in mainlb to be zeroed */
	wa = (double *) malloc((2*m*n+4*n + 11*m*m+8*m) * sizeof(double));
	iwa = (int *) malloc(3*n * sizeof(int));
	strcpy(task, "START");

	while(1){
		/* Main workhorse setulb() from ../appl/lbfgsb.c : */
		setulb(n, m, x, l, u, nbd, &f, g, factr, &pgtol, wa, iwa, task, tr, lsave, isave, dsave);
		/*    Rprintf("in lbfgsb - %s\n", task);*/
		if(strncmp(task, "FG", 2) == 0){
            f = optimGradient(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, n, x, g);
			if(!isfinite(f)){
				cerr << "L-BFGS-B needs finite values of 'fn'" << endl;
				exit(1);
			}
		}else if(strncmp(task, "NEW_X", 5) == 0){
			iter++;
			if(trace == 1 &&(iter % nREPORT == 0)){
				cout << "iter " << iter << " value " << f << endl;
			}
			if(iter > maxit){
				*fail = 1;
				break;
			}
		}else if(strncmp(task, "WARN", 4) == 0){
			*fail = 51;
			break;
		}else if(strncmp(task, "CONV", 4) == 0){
			break;
		}else if(strncmp(task, "ERROR", 5) == 0){
			*fail = 52;
			break;
		}else{ /* some other condition that is not supposed to happen */
			*fail = 52;
			break;
		}
	}

	*Fmin = f;
	*fncount = *grcount = isave[33];
	if(trace){
		cout << "final value " << *Fmin << endl;
		if(iter < maxit && *fail == 0)
			cout << "converged" << endl;
		else
			cout << "stopped after " << iter << " iterations\n";
	}
	strcpy(msg, task);

	free(g);
	free(wa);
	free(iwa);
}


/**
 Function to access the L-BFGS-B function, taken from IQ-TREE package which is further taken from HAL_HAS software package
 1. int nvar : The number of the variables
 2. double* vars : initial values of the variables
 3. double* lower : lower bounds of the variables
 4. double* upper : upper bounds of the variables
 5. double pgtol: gradient tolerance
 5. int maxit : max # of iterations
 @return minimized function value
 After the function is invoked, the values of x will be updated
*/
double L_BFGS_B(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, double* x, double* l, double* u){
  int debug = 0;

  int i;
  double Fmin;
  int fail;
  int fncount;
  int grcount;
  char msg[100];

  int m = 10;          // number of BFGS updates retained in the "L-BFGS-B" method. It defaults to 5.

  int *nbd;           // 0: unbounded; 1: lower bounded; 2: both lower & upper; 3: upper bounded
  nbd = new int[n];
  for(i = 0; i < n; i++) nbd[i] = 2;

  double factr = 1e+7; // control the convergence of the "L-BFGS-B" method.
  // Convergence occurs when the reduction in the object is within this factor
  // of the machine tolerance.
  // Default is 1e7, that is a tolerance of about 1e-8

  double pgtol = opt_type.tolerance;
  double maxit = opt_type.miter;

  //	double pgtol = 0;   // helps control the convergence of the "L-BFGS-B" method.
  //    pgtol = 0.0;
  // It is a tolerance on the projected gradient in the current search direction.
  // Default is zero, when the check is suppressed

  int trace = 0;      // non-negative integer.
  if(debug)
      trace = 1;
  // If positive, tracing information on the progress of the optimization is produced.
  // Higher values may produce more tracing information.

  int nREPORT = 10;   // The frequency of reports for the "L-BFGS-B" methods if "trace" is positive.
  // Defaults to every 10 iterations.

  /*#ifdef USE_OLD_PARAM
  lbfgsb(n, m, x, l, u, nbd, &Fmin, fn, gr1, &fail, ex,
  	factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
  #else*/

  // cout << "initial values in bfgs: ";
  // for(int mi = 0; mi < 3; mi++){
  //     cout << "\t" << x[mi];
  // }
  // cout << "\n";
  //
  // cout << "upper bound in bfgs: ";
  // for(int mi = 0; mi < 3; mi++){
  //     cout << "\t" << u[mi];
  // }
  // cout << "\n";

  lbfgsb(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, n, m, x, l, u, nbd, &Fmin, &fail, factr, pgtol, &fncount, &grcount, maxit, msg, trace, nREPORT);
  //#endif

  if(fail == 51 || fail == 52){
      cout << msg << endl;
  }

  delete[] nbd;

  return Fmin;
}


void max_likelihood_BFGS(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, 
    LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double &min_nlnl, int debug){
    int model = lnl_type.model;
    int cn_max = lnl_type.cn_max;
    int only_seg = lnl_type.only_seg;
    int is_total = lnl_type.is_total;
    int patient_age = lnl_type.patient_age;

    int cons = lnl_type.cons;
    int maxj = opt_type.maxj;
    int opt_one_branch = opt_type.opt_one_branch;

    // Set variables
    int npar_ne = 0;    // number of parameters to estimate
    if(opt_one_branch){
      if(debug){
        cout << "update only one branch " << rtree.current_eid + 1 << endl;
      }
      npar_ne = 1;
    }else{
      if(!cons){  // only estimate internal branches
          npar_ne = 2 * rtree.nleaf - 3;
      }else{    // estimate internal branches and mutation rate
          npar_ne = rtree.nleaf - 1;
      }
    }

    if(debug){
      cout << "\nThere are " << npar_ne << " parameters to estimate " << endl;
    }

    int ndim = get_ndim(maxj, npar_ne, model, only_seg);
    double* variables = new double[ndim + 1];
    memset(variables, 0.0, (ndim + 1) * sizeof(double));
    double* upper_bound = new double[ndim + 1];
    memset(upper_bound, 0.0, (ndim + 1) * sizeof(double));
    double* lower_bound = new double[ndim + 1];
    memset(lower_bound, 0.0, (ndim + 1) * sizeof(double));

    if(cons){    // edges converted to ratio to incorporate time constraints
        // check tip validity
        if(!is_tip_age_valid(rtree.get_node_ages(), opt_type.tobs)){
          cout << "Tip timings inconsistent with observed data when doing BFGS!" << endl;
          rtree.print();
          cout << rtree.make_newick() << endl;         
          exit(1);
        }

        vector<double> ratios = rtree.get_ratio_from_age();

        if(debug){
            cout << "\nInitializing variables related to branch length" << endl;
            cout << "time ratios obtained from node ages: ";
            for(int i = 0; i < rtree.nleaf - 1; i++){
                cout << "\t" << ratios[i];
            }
            cout << endl;
        }

        if(opt_one_branch){   // only optimize one branch
            int nid = rtree.edges[rtree.current_eid].end - rtree.root_node_id;
            int idx = 1;
            variables[idx] = ratios[nid];
            lower_bound[idx] = MIN_RATIO;
            upper_bound[idx] = MAX_RATIO;
        }else{
            // age of root
            variables[1] = ratios[0];
            // make minimal age of root much larger than maximum observed time to avoid very small branch lengths and failing in optimization
            lower_bound[1] = lnl_type.max_tobs * opt_type.scale_tobs + rtree.nleaf * BLEN_MIN;
            // age at 1st sample, so need to add time until last sample
            upper_bound[1] = lnl_type.max_tobs + patient_age;

            for(int i = 1; i < npar_ne; ++i){
              variables[i + 1] = ratios[i];
              lower_bound[i + 1] = MIN_RATIO;
              upper_bound[i + 1] = MAX_RATIO;
            }
        }
    }else{
        if(opt_one_branch){
            variables[1] = rtree.edges[rtree.current_eid].length;
            lower_bound[1] = BLEN_MIN;
            upper_bound[1] = patient_age;
        }else{
            for(int i = 0; i < npar_ne; ++i){
              variables[i + 1] = rtree.edges[i].length;
              lower_bound[i + 1] = BLEN_MIN;
              upper_bound[i + 1] = patient_age;
            }
        }
    }

    if(maxj){    // estimate mutation rates
        if(model == MK){
            int i = npar_ne;
            variables[i + 1] = rtree.mu;
            lower_bound[i + 1] = MIN_MRATE;
            upper_bound[i + 1] = MAX_MRATE;
        }else{
            int i = npar_ne;
            variables[i + 1] = rtree.dup_rate;
            lower_bound[i + 1] = MIN_MRATE;
            upper_bound[i + 1] = MAX_MRATE;

            i = npar_ne + 1;
            variables[i + 1] = rtree.del_rate;
            lower_bound[i + 1] = MIN_MRATE;
            upper_bound[i + 1] = MAX_MRATE;

            if(!only_seg){
                i = npar_ne + 2;
                variables[i + 1] = rtree.chr_gain_rate;
                lower_bound[i + 1] = MIN_MRATE;
                upper_bound[i + 1] = MAX_MRATE;

                i = npar_ne + 3;
                variables[i + 1] = rtree.chr_loss_rate;
                lower_bound[i + 1] = MIN_MRATE;
                upper_bound[i + 1] = MAX_MRATE;

                i = npar_ne + 4;
                variables[i + 1] = rtree.wgd_rate;
                lower_bound[i + 1] = MIN_MRATE;
                upper_bound[i + 1] = MAX_MRATE;
            }
        }
    }

    if(debug){
        if(opt_one_branch){
            cout << "only optimize one branch" << endl;
        }
        cout << "Variables to estimate: " << endl;
        for(int i = 0; i < ndim + 1; i++){
            cout << "\t" << variables[i] << "\tLB:" << lower_bound[i] << "\tUB:" << upper_bound[i] << endl;
        }
    }

    // variables contains the parameters to estimate(branch length, mutation rate)
    min_nlnl = L_BFGS_B(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, ndim, variables + 1, lower_bound + 1, upper_bound + 1);

    // Check the validity of the tree
    if(cons && !is_tree_valid(rtree, lnl_type.max_tobs, patient_age, cons)){
      cout << "The optimized tree after BFGS " << rtree.make_newick() << " is not valid!" << endl;
      exit(1);
    }

    if(debug){
      cout << "log likelihood of current ML tree " << rtree.make_newick() << " is " << -min_nlnl << endl;
    }

    delete [] lower_bound;
    delete [] upper_bound;
    delete [] variables;
}



// Optimizing the branch length of(node1, node2) with BFGS to incorporate constraints imposed by patient age and tip timings.
// When any branch length is updated, the neighbour length needs to be updated
// Optimize mutation rates if necessary
double optimize_one_branch_BFGS(evo_tree& rtree, map<int, vector<vector<int>>>& vobs, OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, Node* node1, Node* node2){
    int debug = 0;
    if(debug){
        cout << "\tOptimizing the branch " << node1->id + 1 << ", " << node2->id + 1 << endl;
    }

    // does not optimize virtual branch from root
    assert(!((node1->id == rtree.nleaf && node2->id == rtree.nleaf - 1) || (node2->id == rtree.nleaf && node1->id == rtree.nleaf - 1)));
        
    rtree.current_it = (Neighbor*) node1->findNeighbor(node2);
    assert(rtree.current_it);
    rtree.current_it_back = (Neighbor*) node2->findNeighbor(node1);
    assert(rtree.current_it_back);
    assert(rtree.current_it->length == rtree.current_it_back->length);

    int eid = rtree.get_edge_id(node1->id, node2->id);
    assert(eid >= 0);
    rtree.current_eid = eid;  // used to track the branch to optimize

    if(debug){
        cout << "tree before optimization by BFGS " << rtree.make_newick() << endl;
        rtree.print();
        cout << "\tUsing BFGS method to optimize the likelihood of one branch length of edge " << eid + 1 << endl;
        // cout << "\taddress of rtree before optimization " << &rtree << endl;
        // cout << "\taddress of node1 before optimization " << &(rtree.node1) << endl;
        // cout << "\taddress of current_it before optimization " << rtree.current_it << endl;
    }

    opt_type.opt_one_branch = 1;
    opt_type.maxj = 0;
    double negative_lh = MAX_NLNL;
    // optimize ratio based on NNI branch(branch length, node times and ages have been updated during optimization)
    max_likelihood_BFGS(rtree, vobs, obs_decomp, comps, lnl_type, opt_type, negative_lh);

    // cout << "\taddress of rtree after optimization " << &rtree << endl;
    // cout << "\taddress of current_it after optimization " << rtree.current_it << endl;

    opt_type.opt_one_branch = 0;
    opt_type.maxj = 1;
    
    if(debug){
        cout << "\tmax logl: " << -negative_lh << " optimized branch length " << rtree.edges[eid].length << endl;
        cout << "tree after optimization by BFGS " << rtree.make_newick() << endl;
        rtree.print();
    }

    return -negative_lh;
}

