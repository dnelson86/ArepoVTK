/*
 * voronoi_3db.cpp (parts of voronoi_3d.c for multithreaded environment)
 * dnelson
 */
 
#include <alloca.h>
#include "voronoi_3db.h"

#ifdef ENABLE_AREPO

const int access_triangles[4][3] = {
  {1, 3, 2},
  {0, 2, 3},
  {0, 3, 1},
  {0, 1, 2}
};

const int edge_start[6] = { 0, 0, 0, 1, 1, 2 };
const int edge_end[6] = { 1, 2, 3, 2, 3, 3 };
const int edge_opposite[6] = { 3, 1, 2, 3, 0, 1 };
const int edge_nexttetra[6] = { 2, 3, 1, 0, 2, 0 };
 
void init_clear_auxmesh(tessellation * T)
{
  //point *p;
  int i, n;

  T->Ndp = 0;
  T->Ndt = 0;
  T->Nvf = 0;

	// copy allocation factors
	/*
	T->Indi.AllocFacNdp = AUXMESH_ALLOC_SIZE / 2;
	T->Indi.AllocFacNdt = AUXMESH_ALLOC_SIZE;
	T->Indi.AllocFacNvf = AUXMESH_ALLOC_SIZE;

	T->MaxNdp = (int)T->Indi.AllocFacNdp;
	T->MaxNdt = (int)T->Indi.AllocFacNdt;
	T->MaxNvf = (int)T->Indi.AllocFacNvf;		

	// allocate
	T->VF = static_cast<face*>
								 (mymalloc_movable(T->VF, "VFaux",T->MaxNvf * sizeof(face)));
	T->DP = static_cast<point*>
								 (mymalloc_movable(T->DP, "DPaux", (T->MaxNdp+5) * sizeof(point)));
	T->DP += 5; // leave first five for bounding tetra + infinity
	T->DT = static_cast<tetra*>
								 (mymalloc_movable(T->DT, "DTaux", T->MaxNdt * sizeof(tetra)));
	
	T->DTC = static_cast<tetra_center*>(
									mymalloc_movable(&T->DTC, "AuxDTC", T->MaxNdt * sizeof(tetra_center)));
	T->DTF = static_cast<char*>(
									mymalloc_movable(&T->DTF, "AuxDTF", T->MaxNdt * sizeof(char)));	
	*/
  /* construct all encompassing huge tetrahedron */

  double box = All.BoxSize;
	double tetra_incircle, tetra_sidelength, tetra_height, tetra_face_height;

  tetra_incircle = 1.5 * box;
  tetra_sidelength = tetra_incircle * sqrt(24);
  tetra_height = sqrt(2.0 / 3) * tetra_sidelength;
  tetra_face_height = sqrt(3.0) / 2.0 * tetra_sidelength;

  //point *DP = T->DP;
  tetra *DT = T->DT;
/*
  // first, let's make the points
  DP[-4].x = 0.5 * tetra_sidelength;
  DP[-4].y = -1.0 / 3 * tetra_face_height;
  DP[-4].z = -0.25 * tetra_height;

  DP[-3].x = 0;
  DP[-3].y = 2.0 / 3 * tetra_face_height;
  DP[-3].z = -0.25 * tetra_height;

  DP[-2].x = -0.5 * tetra_sidelength;
  DP[-2].y = -1.0 / 3 * tetra_face_height;
  DP[-2].z = -0.25 * tetra_height;

  DP[-1].x = 0;
  DP[-1].y = 0;
  DP[-1].z = 0.75 * tetra_height;

  for(i = -4; i <= -1; i++)
    {
      DP[i].x += 0.5 * box;
      DP[i].y += 0.5 * box;
      DP[i].z += 0.5 * box;
    }

  for(i = -4, p = &DP[-4]; i < 0; i++, p++)
    {
      p->index = -1;
      p->task = ThisTask;
      p->timebin = 0;
    }

  // we also define a neutral element at infinity
  DP[DPinfinity].x = GSL_POSINF;
  DP[DPinfinity].y = GSL_POSINF;
  DP[DPinfinity].z = GSL_POSINF;
  DP[DPinfinity].index = -1;
  DP[DPinfinity].task = ThisTask;
  DP[DPinfinity].timebin = 0;
	*/
  // now let's make the big tetrahedron
  DT[0].p[0] = -4;  DT[0].p[1] = -3;  DT[0].p[2] = -2;  DT[0].p[3] = -1;

  for(i = 0; i < 4; i++)
    {
      n = i + 1;		// tetra index

      DT[0].t[i] = n;
      DT[0].s[i] = 3;

      DT[n].t[3] = 0;
      DT[n].s[3] = i;
      DT[n].p[3] = DPinfinity;
    }

  DT[1].p[0] = DT[0].p[1];  DT[1].p[1] = DT[0].p[2];  DT[1].p[2] = DT[0].p[3];
  DT[2].p[0] = DT[0].p[0];  DT[2].p[1] = DT[0].p[3];  DT[2].p[2] = DT[0].p[2];
  DT[3].p[0] = DT[0].p[0];  DT[3].p[1] = DT[0].p[1];  DT[3].p[2] = DT[0].p[3];
  DT[4].p[0] = DT[0].p[0];  DT[4].p[1] = DT[0].p[2];  DT[4].p[2] = DT[0].p[1];

  DT[1].t[0] = 2;  DT[2].t[0] = 1;  DT[1].s[0] = 0;  DT[2].s[0] = 0;
  DT[1].t[1] = 3;  DT[3].t[0] = 1;  DT[1].s[1] = 0;  DT[3].s[0] = 1;
  DT[1].t[2] = 4;  DT[4].t[0] = 1;  DT[1].s[2] = 0;  DT[4].s[0] = 2;
  DT[2].t[2] = 3;  DT[3].t[1] = 2;  DT[2].s[2] = 1;  DT[3].s[1] = 2;
  DT[2].t[1] = 4;  DT[4].t[2] = 2;  DT[2].s[1] = 2;  DT[4].s[2] = 1;
  DT[3].t[2] = 4;  DT[4].t[1] = 3;  DT[3].s[2] = 1;  DT[4].s[1] = 2;

  T->Ndt = 5;			/* we'll start out with 5 tetras */

  CentralOffsetX = 0.5 * box - 0.5000001 * tetra_sidelength;
  CentralOffsetY = 0.5 * box - (1.0000001 / 3) * tetra_face_height;
  CentralOffsetZ = 0.5 * box - 0.25000001 * tetra_height;

  ConversionFac = 1.0 / (1.001 * tetra_sidelength);

  for(i = -4; i < 0; i++)
    set_integers_for_point(T, i);
}

int insert_point_new(tessellation * T, int pp, int ttstart)	/* returns a tetra that (currently) contains the point pp */
{
  int tt0, tt1, tt2, tt3, tt4, tetra_with_p, tt;
  int to_check[STACKSIZE_TETRA], freestack[STACKSIZE_TETRA];
  int n_faces_to_check = 0, nfree_on_stack = 0, moves;
  int tip_index, flag, edgeface_nr;
  int non_convex, convex_edge = 0, i, j;

  /* first, need to do a point location */
  tt0 = get_tetra(T, &T->DP[pp], &moves, ttstart, &flag, &edgeface_nr);

  tetra_with_p = tt0;

  if(flag == 1)			/* that's the normal split of a tetrahedron into 4 */
    {
      if(n_faces_to_check >= STACKSIZE_TETRA - 4)
	terminate("stacksize exceeded");

      /* we now need to split this tetrahedron into four  */
      if(nfree_on_stack)
	tt1 = freestack[--nfree_on_stack];
      else
	tt1 = T->Ndt++;

      if(nfree_on_stack)
	tt2 = freestack[--nfree_on_stack];
      else
	tt2 = T->Ndt++;

      if(nfree_on_stack)
	tt3 = freestack[--nfree_on_stack];
      else
	tt3 = T->Ndt++;

      if(T->Ndt > T->MaxNdt)
				terminate("Ndt > MaxNdt");

      make_a_1_to_4_flip(T, pp, tt0, tt1, tt2, tt3);

      /* now we have a triangulation again - need to check whether there are
         facets that are not Delaunay */

      /* let's initialize a stack with the facets that we need to check */

      n_faces_to_check = 0;

      to_check[n_faces_to_check++] = tt0;
      to_check[n_faces_to_check++] = tt1;
      to_check[n_faces_to_check++] = tt2;
      to_check[n_faces_to_check++] = tt3;
      char *DTF = T->DTF;
      DTF[tt0] = 0;
      DTF[tt1] = 0;
      DTF[tt2] = 0;
      DTF[tt3] = 0;
    }

  if(flag == 2)
    {
      /* create four new tetra  */
      if(nfree_on_stack)
	tt1 = freestack[--nfree_on_stack];
      else
	tt1 = T->Ndt++;

      if(nfree_on_stack)
	tt2 = freestack[--nfree_on_stack];
      else
	tt2 = T->Ndt++;

      if(nfree_on_stack)
	tt3 = freestack[--nfree_on_stack];
      else
	tt3 = T->Ndt++;

      if(nfree_on_stack)
	tt4 = freestack[--nfree_on_stack];
      else
	tt4 = T->Ndt++;

      if(T->Ndt > T->MaxNdt)
				terminate("Ndt > MaxNdt");

      n_faces_to_check = 0;

      to_check[n_faces_to_check++] = tt0;
      to_check[n_faces_to_check++] = T->DT[tt0].t[edgeface_nr];
      to_check[n_faces_to_check++] = tt1;
      to_check[n_faces_to_check++] = tt2;
      to_check[n_faces_to_check++] = tt3;
      to_check[n_faces_to_check++] = tt4;

      char *DTF = T->DTF;
      DTF[tt0] = 0;
      DTF[T->DT[tt0].t[edgeface_nr]] = 0;
      DTF[tt1] = 0;
      DTF[tt2] = 0;
      DTF[tt3] = 0;
      DTF[tt4] = 0;

      make_a_face_split(T, tt0, edgeface_nr, pp, tt1, tt2, tt3, tt4);
    }

  if(flag == 3)			/* here we need to split an edge */
    {
      int i, j, k, l, ii, jj, kk, ll, m, count;
      int prev, next;

      /* count how many triangles share the edge */
      i = edge_start[edgeface_nr];
      j = edge_end[edgeface_nr];
      k = edge_opposite[edgeface_nr];
      l = edge_nexttetra[edgeface_nr];

      count = 0;
      n_faces_to_check = 0;

      prev = tt0;
      do
	{
	  to_check[n_faces_to_check++] = prev;
	  T->DTF[prev] = 0;

	  tetra *DT = T->DT;
	  next = DT[prev].t[l];

	  for(m = 0, ll = ii = jj = -1; m < 4; m++)
	    {
	      if(DT[next].p[m] == DT[prev].p[k])
		ll = m;
	      if(DT[next].p[m] == DT[prev].p[i])
		ii = m;
	      if(DT[next].p[m] == DT[prev].p[j])
		jj = m;
	    }

	  if(ll < 0 || ii < 0 || jj < 0)
	    terminate("inconsistency");

	  kk = 6 - (ll + ii + jj);

	  prev = next;
	  i = ii;
	  l = ll;
	  j = jj;
	  k = kk;

	  count++;

	  if(count > 1000)
	    terminate("count exceeded");
	}
      while(next != tt0);

      //int *ttlist = mymalloc_movable(&ttlist, "ttlist", count * sizeof(int));
			int *ttlist = (int *)alloca(sizeof(int) * count);
			//int ttl[100]; int *ttlist = &ttl[0];

      for(i = 0; i < count; i++)
	{
	  if(nfree_on_stack)
	    ttlist[i] = freestack[--nfree_on_stack];
	  else
	    {
	      ttlist[i] = T->Ndt++;

	      if(T->Ndt > T->MaxNdt)
					terminate("Ndt > MaxNdt");
	    }

	  to_check[n_faces_to_check++] = ttlist[i];
	  T->DTF[ttlist[i]] = 0;
	}
      make_an_edge_split_new(T, tt0, edgeface_nr, count, pp, ttlist);

      //myfree(ttlist);
    }

  int iter = 0;

  while(n_faces_to_check)
    {
      iter++;
      if(iter > 200000)
	terminate("too many iterations");

      tt = to_check[--n_faces_to_check];	/* this is the current tetra to look at.
						   The facet in question lies opposite to q */
      if(T->DT[tt].t[0] < 0)	/* deleted? */
	continue;

      for(tip_index = 0; tip_index < 4; tip_index++)
	if(T->DT[tt].p[tip_index] == pp)
	  break;

      if(tip_index < 4)		/* otherwise the facet has been removed in a 3-2 flip */
	{
	  tetra *DT = T->DT;
	  point *DP = T->DP;
	  int qq = DT[tt].t[tip_index];	/* tetrahedron that's opposite of ours and shares the facet */
	  int ppp = DT[qq].p[DT[tt].s[tip_index]];	/* point that's opposite of the facet in the other tetrahedron */

	  int ret, ret_exact;

	  ret = InSphere_Errorbound(&DP[DT[qq].p[0]], &DP[DT[qq].p[1]], &DP[DT[qq].p[2]], &DP[DT[qq].p[3]], &DP[pp]);
	  //CountInSphereTests++;

	  if(ret != 0)
	    ret_exact = ret;
	  else
	    {
	      // let's decide with exact integer arithmetic
	      ret_exact = InSphere_Exact(&DP[DT[qq].p[0]], &DP[DT[qq].p[1]], &DP[DT[qq].p[2]], &DP[DT[qq].p[3]], &DP[pp]);
	      //CountInSphereTestsExact++;
	    }

	  if(ret_exact > 0)	/* facet is illegal, because point lies inside */
	    {
	      /* let's see whether the point lies in the triangle, or on a side, or opposite of one convex edge */

	      non_convex = convex_edge_test(T, tt, tip_index, &convex_edge);

	      if(non_convex == 0)	/* we can make a 2-3 flip */
		{
		  int ww;

		  if(nfree_on_stack)
		    ww = freestack[--nfree_on_stack];
		  else
		    ww = T->Ndt++;

		  if(T->Ndt > T->MaxNdt)
				terminate("Ndt > MaxNdt");

		  if(n_faces_to_check >= STACKSIZE_TETRA - 3)
		    terminate("stacksize exceeded");

		  make_a_2_to_3_flip(T, tt, tip_index, qq, T->DT[tt].s[tip_index], ppp, ww);

		  to_check[n_faces_to_check++] = tt;
		  to_check[n_faces_to_check++] = qq;
		  to_check[n_faces_to_check++] = ww;
		  T->DTF[tt] = 0;
		  T->DTF[qq] = 0;
		  T->DTF[ww] = 0;
		}
	      else if(non_convex == 1)	/* we might be able to make a 3-2 flip, or we deal with a convex edge on the outer hull */
		{
		  /* test whether the reflex edge is surrounded by exactly three tetrahedra */

		  i = convex_edge + 2;
		  if(i >= 3)
		    i -= 3;
		  i = access_triangles[tip_index][i];

		  for(j = 0; j < 4; j++)
		    if(DT[tt].p[i] == DT[qq].p[j])
		      break;

		  if(j >= 4)
		    {
		      terminate("not found");
		    }


		  if(DT[tt].t[i] == DT[qq].t[j])	/* this means there is exactly one tetrahedron between them, i.e. we have found the
							   third partner for the flip */
		    {
		      int ww;

		      ww = DT[tt].t[i];

		      make_a_3_to_2_flip(T, tt, qq, ww, tip_index, convex_edge, DT[tt].s[tip_index]);

		      DT[ww].t[0] = -1;	/* mark as deleted */

		      if(nfree_on_stack < STACKSIZE_TETRA)
			freestack[nfree_on_stack++] = ww;
		      else
			terminate("stack full");


		      tetra_with_p = tt;
		      if(n_faces_to_check >= STACKSIZE_TETRA - 2)
			terminate("stack too full");

		      to_check[n_faces_to_check++] = tt;
		      to_check[n_faces_to_check++] = qq;
		      T->DTF[tt] = 0;
		      T->DTF[qq] = 0;
		    }
		  else
		    {
		      if(DT[DT[tt].t[i]].p[DT[tt].s[i]] == DPinfinity && DT[DT[qq].t[j]].p[DT[qq].s[j]] == DPinfinity)
			{
			  printf("convex edge between points=%d %d on outer hull found\n",
				 (int) (DT[tt].p[access_triangles[tip_index][convex_edge]]),
				 (int) (DT[tt].p[access_triangles[tip_index][convex_edge < 2 ? convex_edge + 1 : 0]]));

			  terminate("inconsistency");	/* this should not occur since we have embedded the points into a convex big triangle */
			}
		    }
		}
	      else if(non_convex == 2)	/* we might be able to make a 4-4 flip */
		{
		  i = convex_edge + 2;
		  if(i >= 3)
		    i -= 3;
		  i = access_triangles[tip_index][i];	/* this is the point opposite of edge (but not tip) */

		  tetra *DT = T->DT;
		  char *DTF = T->DTF;

		  for(j = 0; j < 4; j++)
		    if(DT[tt].p[i] == DT[qq].p[j])
		      break;

		  if(DT[DT[tt].t[i]].p[DT[tt].s[i]] == DT[DT[qq].t[j]].p[DT[qq].s[j]])
		    {
		      /* ok, so we really have 4 tetra. The opposite points match up */
		      //                      printf("we should do a 4-4 flip\n");

		      to_check[n_faces_to_check++] = tt;
		      to_check[n_faces_to_check++] = qq;
		      to_check[n_faces_to_check++] = DT[tt].t[i];
		      to_check[n_faces_to_check++] = DT[qq].t[j];
		      DTF[tt] = 0;
		      DTF[qq] = 0;
		      DTF[DT[tt].t[i]] = 0;
		      DTF[DT[qq].t[j]] = 0;

		      make_a_4_to_4_flip(T, tt, tip_index, convex_edge);
		    }
		}
	    }
	  else
	    tetra_with_p = tt;
	}
    }

  return tetra_with_p;
}

void make_an_edge_split_new(tessellation * T, int tt0, int edge_nr, int count, int pp, int *ttlist)
{
  tetra *DT = T->DT;
  tetra *t0 = &DT[tt0];
  tetra *prev, *next;
  tetra **tlist, **t_orig_list;
  int *i_list, *j_list, *k_list, *l_list;
  int i, j, k, l, ii, jj, kk, ll, m, nr, nrm, nrp;

  //Count_EdgeSplits++;
  //CountFlips++;
/*
  tlist = mymalloc("tlist", count * sizeof(tetra *));
  t_orig_list = mymalloc("t_orig_list", count * sizeof(tetra *));
  i_list = mymalloc("i_list", sizeof(int) * count);
  j_list = mymalloc("j_list", sizeof(int) * count);
  k_list = mymalloc("k_list", sizeof(int) * count);
  l_list = mymalloc("l_list", sizeof(int) * count);
*/	

	tlist = (tetra **)alloca(sizeof(tetra *) * count);
	t_orig_list = (tetra **)alloca(sizeof(tetra *) * count);
	
	i_list = (int *)alloca(sizeof(int) * count);
	j_list = (int *)alloca(sizeof(int) * count);
	k_list = (int *)alloca(sizeof(int) * count);
	l_list = (int *)alloca(sizeof(int) * count);

  for(i = 0; i < count; i++)
    tlist[i] = &DT[ttlist[i]];

  i = edge_start[edge_nr];
  j = edge_end[edge_nr];
  k = edge_opposite[edge_nr];
  l = edge_nexttetra[edge_nr];

  nr = 0;
  prev = t0;
  do
    {
      t_orig_list[nr] = prev;
      i_list[nr] = i;
      j_list[nr] = j;
      k_list[nr] = k;
      l_list[nr] = l;

      next = &DT[prev->t[l]];

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
	{
	  if(next->p[m] == prev->p[k])
	    ll = m;
	  if(next->p[m] == prev->p[i])
	    ii = m;
	  if(next->p[m] == prev->p[j])
	    jj = m;
	}

      if(ll < 0 || ii < 0 || jj < 0)
	terminate("inconsistency");

      kk = 6 - (ll + ii + jj);

      prev = next;
      i = ii;
      l = ll;
      j = jj;
      k = kk;

      nr++;
    }
  while(next != t0);


  for(nr = 0; nr < count; nr++)
    {
      *tlist[nr] = *t_orig_list[nr];

      t_orig_list[nr]->p[j_list[nr]] = pp;
      tlist[nr]->p[i_list[nr]] = pp;


      t_orig_list[nr]->t[i_list[nr]] = tlist[nr] - DT;
      tlist[nr]->t[j_list[nr]] = t_orig_list[nr] - DT;

      t_orig_list[nr]->s[i_list[nr]] = j_list[nr];
      tlist[nr]->s[j_list[nr]] = i_list[nr];

      DT[tlist[nr]->t[i_list[nr]]].t[tlist[nr]->s[i_list[nr]]] = tlist[nr] - DT;

      nrp = nr + 1;
      if(nrp >= count)
	nrp -= count;

      nrm = nr - 1;
      if(nrm < 0)
	nrm += count;

      tlist[nr]->t[l_list[nr]] = tlist[nrp] - DT;
      tlist[nr]->s[l_list[nr]] = k_list[nrp];

      tlist[nr]->t[k_list[nr]] = tlist[nrm] - DT;
      tlist[nr]->s[k_list[nr]] = l_list[nrm];
    }

		/*
  myfree(l_list);
  myfree(k_list);
  myfree(j_list);
  myfree(i_list);
  myfree(t_orig_list);
  myfree(tlist); */
}

void compute_auxmesh_volumes(tessellation *T, double *vol)
{
  int i, bit, nr;

  for(i = 0; i < T->Ndp; i++)
    vol[i] = 0;

  //Edge_visited = mymalloc_movable(&Edge_visited, "Edge_visited", T->Ndt * sizeof(unsigned char)); // Edge_visited global
	unsigned char *visited_edges = (unsigned char *)alloca(sizeof(unsigned char) * T->Ndt);

  for(i = 0; i < T->Ndt; i++)
    visited_edges[i] = 0;

  for(i = 0; i < T->Ndt; i++)
    {
      if(T->DT[i].t[0] < 0)      /* deleted ? */
        continue;

      bit = 1;
      nr = 0;

      while(visited_edges[i] != EDGE_ALL)
        {
          if((visited_edges[i] & bit) == 0)
            derefine_refine_process_edge_new(T, vol, i, nr, visited_edges);

          bit <<= 1;
          nr++;
        }
    }

  //myfree(Edge_visited);
}

void derefine_refine_process_edge_new(tessellation * T, double *vol, int tt, int nr, unsigned char *visited_edges)
{
  tetra *DT = T->DT;
  point *DP = T->DP;
  tetra_center *DTC = T->DTC;

  int i, j, k, l, m, ii, jj, kk, ll, nn, count, nr_next, p1, p2;
  tetra *prev, *next;
  tetra_center *prevc, *nextc;
  double ax, ay, az;
  double bx, by, bz;
  double cx, cy, cz;
  double nx, ny, nz;
  double hhx, hhy, hhz;
  double darea, dvol, h;

  tetra *t = &DT[tt];

  i = edge_start[nr];
  j = edge_end[nr];
  k = edge_opposite[nr];
  l = edge_nexttetra[nr];

  visited_edges[tt] |= (1 << nr);

  p1 = t->p[i];
  p2 = t->p[j];

  double area = 0;

  cx = DTC[tt].cx;
  cy = DTC[tt].cy;
  cz = DTC[tt].cz;

  count = 0;

  prev = t;
  prevc = &DTC[tt];
  do
    {
      nn = prev->t[l];
      next = &DT[nn];
      nextc = &DTC[nn];

      if(prev != t && next != t)
	{
	  ax = prevc->cx - cx;
	  ay = prevc->cy - cy;
	  az = prevc->cz - cz;

	  bx = nextc->cx - cx;
	  by = nextc->cy - cy;
	  bz = nextc->cz - cz;

	  nx = ay * bz - az * by;
	  ny = az * bx - ax * bz;
	  nz = ax * by - ay * bx;

	  darea = 0.5 * sqrt(nx * nx + ny * ny + nz * nz);
	  area += darea;
	}

      for(m = 0, ll = ii = jj = -1; m < 4; m++)
	{
	  if(next->p[m] == prev->p[k])
	    ll = m;
	  if(next->p[m] == prev->p[i])
	    ii = m;
	  if(next->p[m] == prev->p[j])
	    jj = m;
	}

      if(ll < 0 || ii < 0 || jj < 0)
	terminate("inconsistency");

      kk = 6 - (ll + ii + jj);


      /* need to determine the edge number to be able to flag it */

      for(nr_next = 0; nr_next < 6; nr_next++)
	if((edge_start[nr_next] == ii && edge_end[nr_next] == jj) || (edge_start[nr_next] == jj && edge_end[nr_next] == ii))
	  {
	    if((visited_edges[nn] & (1 << nr_next)) && next != t)
	      terminate("inconsistency");

	    visited_edges[nn] |= (1 << nr_next);
	    break;
	  }

      prev = next;
      prevc = nextc;
      i = ii;
      l = ll;
      j = jj;
      k = kk;

      count++;

      if(count > 1000)
	terminate("count is too large");
    }
  while(next != t);

  i = edge_start[nr];
  j = edge_end[nr];

  hhx = 0.5 * (DP[p1].x - DP[p2].x);
  hhy = 0.5 * (DP[p1].y - DP[p2].y);
  hhz = 0.5 * (DP[p1].z - DP[p2].z);

  h = sqrt(hhx * hhx + hhy * hhy + hhz * hhz);
  dvol = (1.0 / 3) * area * h;

  if(p1 >= 0 && p1 < T->Ndp)
    vol[p1] += dvol;

  if(p2 >= 0 && p2 < T->Ndp)
    vol[p2] += dvol;
}

#endif // ENABLE_AREPO
