from numpy import *
import copy
from astropy.io import fits
import time
from scipy.sparse import lil_matrix

# Version 1: Starting to keep track of versions
# 1.01: added wmat_to_rho
# 1.02: Prints out diagnostics for incorrectly-sized Wmat in LM
# 1.03: Added 'verbose' option to secderiv
# 1.04: Added wmat_to_errormat
# 1.05: Added option to save jacobian from LM as something
# 1.06: Added wrapper around miniLM, including cov mat return
# 1.07: Added wrapper around miniNM
# 1.08: Added option to turn Cmat computation off for miniNM_new
# 1.10: Jacobian now stored as sparse matrix.
# 1.11: 'verbose' option for secderiv passed through miniNM_new
# 1.12: Added 'save_patches'
# 1.13: Added nan check to miniLM
# 1.14: Added option to miniLM to use dense jacobian matrices (faster, if they are dense)
# 1.15: Added option to miniLM_new to not return Cmat, won't recompute Jacobian after search!
# 1.16: Added option to miniLM to use multiprocessing pool
# 1.17: Added eigenvector decomposition (better linalg.eig)
# 1.18: Added save_jacobian to miniLM
# 1.19: pool works now in miniLM
# 1.20: Warns about determinant
# 1.21: Secder automatic sign flip if it sees the fit is up against a bound
# 1.22: Improved version of save_patches

version = 1.22

print(f"DavidsNM Version {version}")

def eig(the_matrix):
    if any(isnan(the_matrix)):
        print("Couldn't decompose matrix!")
        return 0, 0, 0

    evals, evecs = linalg.eig(the_matrix)

    inds = argsort(evals)[::-1]

    evecs = transpose(evecs)
    evecs = evecs[inds]
    evals = evals[inds]
    evecs_norm = array([evecs[i]*sqrt(evals[i]) for i in range(len(evecs))])
    return evals, evecs, evecs_norm


def save_img(dat, imname):

    fitsobj = fits.HDUList()
    hdu = fits.PrimaryHDU()
    hdu.data = dat
    fitsobj.append(hdu)
    fitsobj.writeto(imname, overwrite=True)
    fitsobj.close()



def save_patches(patch_list, imname):
    try:
        patch1 = len(patch_list[0])
        patch2 = len(patch_list[0][0])
    except:
        save_img(zeros([0,0]), imname)
        return 0

    ideal_pix_on_a_side = sqrt(len(patch_list)*patch1*patch2*1.)
    grid1 = int(ceil(ideal_pix_on_a_side/patch1))
    grid2 = int(ceil(ideal_pix_on_a_side/patch2))

    all_data = zeros([grid1*patch1, grid2*patch2], dtype=float64) + sqrt(-1.)

    max_i = 0
    max_j = 0
    for i in range(len(patch_list)):
        ipos = i % grid1
        jpos = int(floor(i / grid1))

        all_data[ipos*patch1 : (ipos + 1)*patch1,
                 jpos*patch2 : (jpos + 1)*patch2] = patch_list[i]

        max_i = max(max_i, (ipos + 1)*patch1)
        max_j = max(max_j, (jpos + 1)*patch2)

    save_img(all_data[:max_i, :max_j], imname)


def get_simpons_weights(n):
    if n%2 != 1:
        print("n should be odd!", n)
        return None
    weights = zeros(n, dtype=float64)
    weights[1::2] = 4.
    weights[0::2] = 2.
    weights[0] = 1.
    weights[n - 1] = 1.

    return weights/sum(weights)


# Start of minimization routines

def wmat_to_rhomat(wmat):
    covmat = linalg.inv(wmat)
    rhomat = zeros(wmat.shape, dtype=float64)


    for i in range(len(wmat)):
        for j in range(len(wmat)):
            rhomat[i,j] = covmat[i,j]/sqrt(covmat[i,i]*covmat[j,j])
    return rhomat

def wmat_to_errormat(wmat):
    covmat = linalg.inv(wmat)
    errormat = zeros(wmat.shape, dtype=float64)


    for i in range(len(wmat)):
        for j in range(len(wmat)):
            errormat[i,j] = sign(covmat[i,j])*sqrt(abs(covmat[i,j]))
    return errormat

def smoothprior(x): # quadratic for small x, linear for large
    return 1000.*(sqrt(1. + x**2.) - 1.)



def f(P, merged_list):
    chi2fn = merged_list[0]
    return chi2fn(P, merged_list[2:])


def listsort(pf):
    [P, F] = pf

    tmplist = []
    for i in range(len(P)):
        tmplist.append([F[i], P[i].tolist()])
    tmplist.sort()

    P = []
    F = []
    for i in range(len(tmplist)):
        P.append(tmplist[i][1])
        F.append(tmplist[i][0])
    P = array(P, dtype=float64)
    return [P, F]

def strlist(list):
    newlist = []
    for i in range(len(list)):
        newlist.append(str(list[i]))
    return newlist

def find_best_dx(displ_params, free_params, goal, i, P, minchi, merged_list, verbose, sign_to_use = 1.,):
    dx_scale_max = -1e6
    dx_scale_min = 1e6

    dx_scale = 1.
    dP = zeros(len(P[0]), dtype=float64)
    triesfordx = 0

    dP[i] = displ_params[free_params.index(i)]*dx_scale*sign_to_use
    fPdP = f(P[0] + dP, merged_list)

    while abs(fPdP - minchi) > 2.*goal or abs(fPdP - minchi) < goal/2.:
        if verbose:
            print(i, dx_scale, dx_scale_min, dx_scale_max, fPdP)

        if abs(fPdP - minchi) < goal/2.: # point is too close to chi2 minimum
            dx_scale_min = min([dx_scale_min, dx_scale])
            dx_scale *= 2.0

        if abs(fPdP - minchi) > 2.*goal: # point is too far away from chi2 minmum
            dx_scale_max = max([dx_scale_max, dx_scale])
            dx_scale /= 2.0

        if dx_scale_min < dx_scale_max:
            dx_scale = (dx_scale_max + dx_scale_min)/2.


        dP[i] = displ_params[free_params.index(i)]*dx_scale*sign_to_use
        fPdP = f(P[0] + dP, merged_list)

        triesfordx += 1

        if triesfordx > 100:
            print("Couldn't get dx, i = ", i)
            return None
    return displ_params[free_params.index(i)]*dx_scale*sign_to_use


def secderiv(P, merged_list, displ_list, goal, verbose = True): # goal is delta chi2 from minimum, sqrt(goal) is about the sigma
    if verbose:
        print("goal ", goal)

    free_params = []
    displ_params = []
    for i in range(len(displ_list)):
        if displ_list[i] != 0:
            free_params.append(i)
            displ_params.append(displ_list[i]/10.)
    if verbose:
        print("free_params ", free_params)
        print("displ_params ", displ_params)

    W = zeros([len(free_params)]*2, dtype=float64)

    minchi = f(P[0], merged_list)
    dx = zeros(len(free_params), dtype=float64)


    for i in free_params:
        tmp_dx = find_best_dx(displ_params = displ_params, free_params = free_params, goal = goal, i = i, P = P, minchi = minchi, merged_list = merged_list, verbose = verbose, sign_to_use = 1.)
        if tmp_dx == None:
            print("Flipping sign!")
            tmp_dx = find_best_dx(displ_params = displ_params, free_params = free_params, goal = goal, i = i, P = P, minchi = minchi, merged_list = merged_list, verbose = verbose, sign_to_use = -1.)
        if tmp_dx == None:
            print("Tried sign flip, gave up!")
            sys.exit(1)

        dx[free_params.index(i)] = tmp_dx
        if verbose:
            print("dx[free_params.index(i)], ", dx[free_params.index(i)])
    if verbose:
        print("dx ", dx)

    F0 = f(P[0], merged_list)

    Pcollection = [P[0]]
    Fcollection = [f(P[0], merged_list)]

    for i in free_params:
        if verbose:
            print(i)
        for j in free_params:

            if j >= i:
                dP1 = zeros(len(P[0]), dtype=float64)
                dP1[i] = dx[free_params.index(i)]
                dP2 = zeros(len(P[0]), dtype=float64)
                dP2[j] = dx[free_params.index(j)]

                F0 = f(P[0] - (dP1 + dP2)*0.5, merged_list)
                F1 = f(P[0] + dP1 - (dP1 + dP2)*0.5, merged_list)
                F12 = f(P[0] + dP1 + dP2 - (dP1 + dP2)*0.5, merged_list)
                F2 = f(P[0] + dP2 - (dP1 + dP2)*0.5, merged_list)

                Pcollection.append(P[0] + dP1)
                Fcollection.append(F1)
                Pcollection.append(P[0] + dP2)
                Fcollection.append(F2)
                Pcollection.append(P[0] + dP1 + dP2)
                Fcollection.append(F12)

                W[free_params.index(i), free_params.index(j)] = (F12 - F1 - F2 + F0)/(dx[free_params.index(i)]*dx[free_params.index(j)])
                W[free_params.index(j), free_params.index(i)] = W[free_params.index(i), free_params.index(j)]


    W /= 2.0

    if verbose:
        print("Weight Matrix ", W)
    return [W, Pcollection, Fcollection]




def fillindx(dx, free_params, displ_list):
    tmp_dx = zeros(  len(displ_list), dtype=float64)


    for i in range(len(free_params)):
        tmp_dx[free_params[i]] = dx[i]
    return tmp_dx



def better_secderiv(P, merged_list, displ_list, goal): # goal is delta chi2 from minimum, sqrt(goal) is about the sigma
    print("goal ", goal)


    free_params = []
    displ_params = []
    for i in range(len(displ_list)):
        if displ_list[i] != 0:
            free_params.append(i)
            displ_params.append(abs(displ_list[i]))
    print("free_params ", free_params)
    print("displ_params ", displ_params)

    W = zeros([len(free_params)]*2, dtype=float64)

    minchi = f(P[0], merged_list)

    dx = []
    for i in range(len(free_params)):
        dx.append(  zeros(len(free_params), dtype=float64)  )
        dx[i][i] = displ_params[i]






    for k in range(2):


        print("dx ")
        print(dx)

        for i in range(len(free_params)): # Normalize the dxs
            scale = 1. # scale factor for dx[i]

            scale_max = 1.e10
            scale_min = 0.0

            triesforscale = 0


            tmp_dx = fillindx(dx[i], free_params, displ_list)

            fPdP = f(P[0] + scale*tmp_dx, merged_list)

            while abs(fPdP - minchi) > 2.*goal or abs(fPdP - minchi) < goal/2.:
                #print i, tmp_dx, tmp_max, tmp_min, f(P[0] + dP, merged_list)
                if abs(fPdP - minchi) < goal/2.: # too close to minimum
                    scale_min = max(scale_min, scale) # must be at least this far away
                    scale *= 2.0

                if abs(fPdP - minchi) > 2.*goal:
                    scale_max = min(scale_max, scale)
                    scale /= 2.0

                if scale_max != 1.e10 and scale_min != 0.0:
                    scale = (scale_max + scale_min)/2.


                fPdP = f(P[0] + tmp_dx*scale, merged_list)

                if abs(fPdP - minchi) <= 2.*goal and abs(fPdP - minchi) >= goal/2.:
                    dx[i] *= scale

                triesforscale += 1

                if triesforscale > 100:
                    print("Couldn't get dx, i = ", i)
                    sys.exit(1)

            print("dx[i], ", dx[i])
        print("dx ", dx)


        for i in range(len(free_params)):
            for j in range(len(free_params)):

                tmp_dx1 = fillindx(dx[i], free_params, displ_list)
                tmp_dx2 = fillindx(dx[j], free_params, displ_list)

                pt0 = P[0]
                pt1 = P[0] + tmp_dx1
                pt2 = P[0] + tmp_dx2
                pt3 = P[0] + tmp_dx1 + tmp_dx2

                W[i, j] = (
                    f(pt0, merged_list) - f(pt1, merged_list) - f(pt2, merged_list) + f(pt3, merged_list)
                    )/(
                    sqrt(dot(tmp_dx1, tmp_dx1))*
                    sqrt(dot(tmp_dx2, tmp_dx2)))

                W[j, i] = W[i, j]
                # I guess this line explicitly enforces symmetry

        print("Weight Matrix iter ", k)
        print(W)

        eig_vec = linalg.eig(W)[1]
        print("eig_vec")
        print(eig_vec)


        dx = eig_vec



    W = dot(transpose(dx), dot(W, linalg.inv(transpose(dx))))
    W /= 2.0
    print("Weight Matrix ")
    print(W)




    return W





def linfit(x1, y1, x2, y2, targ, limit1, limit2):
    print("linfit ", [x1, y1], [x2, y2])
    slope = (y1 - y2)/(x1 - x2)
    inter = (x1*y2 - x2*y1)/(x1 - x2)


    bestguess = (targ - inter)/slope

    if bestguess > limit1 and bestguess > limit2: # too high
        bestguess = max(limit1, limit2)

    if bestguess < limit1 and bestguess < limit2: # too low
        bestguess = min(limit1, limit2)

    return bestguess

def minos_f(ministarts, minioffsets, merged_list, dx, pos):

    bestF = -1

    tmpstarts = array(ministarts, dtype=float64)

    try:
        len(pos)
        tmpstarts += pos*dx

        if any(tmpstarts*pos != tmpstarts[0]*pos):
            print("tmpstarts are conflicting!")
            print("tmpstarts ", tmpstarts)
            print("ministarts ", ministarts)
            sys.exit(1)
    except:
        tmpstarts[:,pos] += dx

        if any(tmpstarts[:,pos] != tmpstarts[0,pos]):
            print("tmpstarts are conflicting!")
            print("tmpstarts ", tmpstarts)
            print("ministarts ", ministarts)
            sys.exit(1)

    tmpstarts = tmpstarts.tolist()

    for i in range(len(ministarts)):
        [P, F] = miniNM(tmpstarts[i], [1.e-6, 1.e-8], merged_list, minioffsets[i], 0)
        #[P, F] = miniNM(tmpstarts[i], [1.e-8, 1.e-10], merged_list, minioffsets[i], 0)


        if F[0] < bestF or bestF == -1:
            bestF = F[0]
            bestP = P[0]
        if isnan(bestF):
            print("Nan!")
            return [1.e20, P[0]]


    return [bestF, bestP]




def minos(Pmins, minioffsets, chi2fn, dx, targetchi2, pos, minichi2, inlimit = lambda x: 1, passdata = None): #Pmins, minioffsets are lists of starting conditions
    """This is the version to use, not the newer one."""
    merged_list = [chi2fn, inlimit]
    merged_list.extend([passdata])

    toohigh = 0.0
    toolow = 0.0
    toohighchi2 = 0.0
    toolowchi2 = 0.0
    print("targetchi2 ", targetchi2)

    [cur_chi, bestP] = minos_f(Pmins, minioffsets, merged_list, dx, pos)
    if cur_chi < minichi2:
        print("You didn't converge the starting fit!")
        print(bestP)
        return [0., bestP]

    if cur_chi == 1000000:
        print("Minos Start Error!")
        sys.exit(1)

    dxscale = 1.

    chi2list = []

    minostries = 0

    while abs(cur_chi - targetchi2) > 0.000001 and abs(toohigh - toolow) > 0.00001*abs(toohigh) or toohigh == 0.0 or toolow == 0.0:

        dxscale *= 1.5
        print("[abs(cur_chi - targetchi2), max(abs(toohigh - toolow))] ", [abs(cur_chi - targetchi2), abs(toohigh - toolow)])

        chi2list.append([abs(cur_chi - targetchi2), dx, cur_chi])
        chi2list.sort() # low to high

        if cur_chi > targetchi2:
            print(cur_chi, " > ", targetchi2)
            toohigh = dx

            if toohigh == 0.0 or toolow == 0.0:
                dx /= dxscale
            toohighchi2 = cur_chi
        else:
            print(cur_chi, " <= ", targetchi2)
            toolow = dx

            if toohigh == 0.0 or toolow == 0.0:
                dx *= dxscale
            toolowchi2 = cur_chi
        if toohigh != 0.0 and toolow != 0.0:

            print("chi2list [abs(cur_chi - targetchi2), dx, cur_chi] ",chi2list)

            if minostries < 15:
                dx = linfit(chi2list[0][1], chi2list[0][2], chi2list[1][1], chi2list[1][2], targetchi2,
                            0.9*toolow + 0.1*toohigh,
                            0.9*toohigh + 0.1*toolow)
                print("dxlinfit ", dx)
            else:
                dx = (toohigh + toolow)/2.
                print("dx the slow way ", dx)


        [cur_chi, bestP] = minos_f(Pmins, minioffsets, merged_list, dx, pos)
        if cur_chi < minichi2:
            print("You didn't converge the starting fit!")
            print(bestP)
            return [0., bestP]
        print("toohigh, toolow, dx, cur_chi ", toohigh, toolow, dx, cur_chi)
        if isinf(dx):
            return [1.e100, None]
        minostries += 1

    print("Finished getting dx")
    return [dx, bestP]



def better_minos(P0, Pstart, minioffset, merged_list, dx, target_chi2, verbose):
    print("target_chi2 ", target_chi2)

    newmerged_list = copy.deepcopy(merged_list)
    minos_params = [target_chi2, dx, P0, merged_list[0]]
    newmerged_list.append(minos_params)


    newmerged_list[0] = better_minos_chi2fn

    [P, F] = miniNM(Pstart, [1.e-6, 0.], newmerged_list, minioffset, verbose)

    return [dot(P[0] - P0, dx)/sqrt(dot(dx, dx)), P[0] - P0, merged_list[0](P[0], merged_list[2:])]

def better_minos_chi2fn(P, merged_list):
    minos_params = merged_list[-1]

    [target_chi2, dx, P0, chi2fn] = minos_params

    chi2 = chi2fn(P, merged_list[:-1])


    gradient = -dot(P - P0, dx)/dot(dx, dx)
    if abs(gradient) > 1.e2:
        return 0.

    # 1.e5 assures chi2 is positive
    return chi2 + smoothprior(chi2 - target_chi2) + gradient + 1.e5



def improve(pf, merged_list):#P is sorted lowest chi2 to highest
    [P, F] = pf


    M = sum(P[:-1], axis = 0)/len(P[:-1])




    W = P[-1]
    fW = F[-1]

    R = M + M - W
    E = R + (R - M)

    if merged_list[1](E):
        fE = f(E, merged_list)
        if fE < fW:
            P[-1] = E
            F[-1] = fE
            return [P, F]

    if merged_list[1](R):
        fR = f(R, merged_list)
        if fR < fW:
            P[-1] = R
            F[-1] = fR
            return [P, F]


    C1 = 0.5*(M + W)
    C2 = 0.5*(M + R)

    if merged_list[1](C1):
        fC1 = f(C1, merged_list)
        if fC1 < fW:
            P[-1] = C1
            F[-1] = fC1
            return [P, F]

    if merged_list[1](C2):
        fC2 = f(C2, merged_list)
        if fC2 < fW:
            P[-1] = C2
            F[-1] = fC2
            return [P, F]

    for i in range(1,len(P)):
        P[i] = 0.5*(P[0] + P[i])
        F[i] = f(P[i], merged_list)

    return [P, F]





def get_start(P0, displ_list, merged_list):
    P = array([P0]*(len(displ_list) - displ_list.count(0.) + 1), dtype=float64)

    F = [f(P[0], merged_list)]

    j = 1
    for i in range(len(displ_list)):
        if displ_list[i] != 0.:

            P[j,i] += displ_list[i]

            if merged_list[1](P[j]) == 1: # merged_list[1] is inlimit
                F.append(f(P[j], merged_list))
            else:
                print("Changing sign! ", displ_list[i])
                P[j,i] -= 2.*displ_list[i]

                if merged_list[1](P[j]) == 1: # merged_list[1] is inlimit
                    print("Change worked!")
                    F.append(f(P[j], merged_list))
                else:
                    print("start out of range!")
                    print(P[j])
                    return [P, F]
            j += 1


    [P, F] = listsort([P, F])

    return [P, F]



def miniNM(P0, e, merged_list, displ_list, verbose, maxruncount = 15, negativewarning = True, maxiter = 100000):
    runcount = 0

    print("maxruncount ", maxruncount)


    old_F = -1.
    F = [-2.]

    while runcount < maxruncount and old_F != F[0]:
        [P, F] = get_start(P0, displ_list, merged_list)
        if len(F) != len(P): # Started against a limit
            print("Returning starting value!")
            return [array([P0], dtype=float64),
                    array([f(array(P0, dtype=float64), merged_list)], dtype=float64)]
        old_F = F[0]

        k = 1

        noimprovement = 0
        noimprove_F = -1.

        if runcount == 0:
            tmpe = e[0]/10.
            tmpe2 = e[1]/10.
        else:
            tmpe = e[0] # Just checking previous result
            tmpe2 = e[1]


        while (F[-1] > F[0] + tmpe and
               k < maxiter and
               noimprovement < 200 and
               max(abs(  (P[0] - P[-1])/max(max(P[0]), 1.e-10)  )) > tmpe2 ) or k < 2:

            last_F = F[0]
            [P, F] = improve([P, F], merged_list) # Run an iteration
            [P, F] = listsort([P, F])


            if last_F == F[0] or F[0] == old_F:
                noimprovement += 1
            else:
                noimprovement = 0
                tmpe = e[0]/10. # If improvement, run extra
                tmpe2 = e[1]/10.


            if verbose == 1:
                print(P[0], F[0], F[-1] - F[0], k, noimprovement)

            if F[0] < -1.e-8 and negativewarning:
                print("F ", F)
                print("P ", P)
                print("Negative Chi2")
                sys.exit(1.)

            k += 1


        if verbose != -1:
            print("iter F[0]", k, F[0])


        P0 = P[0].tolist()

        runcount += 1

    return [P, F]

def miniNM_new(ministart, miniscale, passdata, chi2fn = None, residfn = None, inlimit = lambda x: True, verbose = False, maxruncount = 15, negativewarning = False, maxiter = 10000, tolerance = [1.e-8, 0], compute_Cmat = True):
    if chi2fn == None:
        chi2fn = lambda x, y: (residfn(x, y)**2.).sum()

    try:
        miniscale = miniscale.tolist()
    except:
        pass

    try:
        ministart = ministart.tolist()
    except:
        pass

    [P, F] = miniNM(ministart, tolerance, [chi2fn, inlimit, passdata], miniscale, verbose = verbose, maxruncount = maxruncount, negativewarning = negativewarning, maxiter = maxiter)

    if compute_Cmat:
        [Wmat, NA, NA] = secderiv(P, [chi2fn, inlimit, passdata], miniscale, 1.e-1, verbose = verbose)
    else:
        Wmat = []

    Cmat = None
    if len(Wmat) > 0:
        if linalg.det(Wmat) != 0:
            Cmat = linalg.inv(Wmat)

    return P[0], F[0], Cmat





def err_from_cov(matrix):
    errs = []
    for i in range(len(matrix)):
        errs.append(sqrt(matrix[i,i]))
    return errs



# Start L-M


def Jacobian(modelfn, unpad_offsetparams, merged_list, unpad_params, displ_list, params, datalen, use_dense_J, pool = None):

    if use_dense_J:
        J = zeros([datalen, len(unpad_params)], dtype=float64, order = 'F')
    else:
        J = zeros([datalen, len(unpad_params)], dtype=float64)#, order = 'F')

    #J = lil_matrix((datalen, len(unpad_params)))

    if pool == None:
        base_mod_list = modelfn(get_pad_params(unpad_params, displ_list, params), merged_list)
    else:
        base_mod_list = modelfn((get_pad_params(unpad_params, displ_list, params), merged_list))

    if pool == None:
        for j in range(len(unpad_params)):
            dparams = copy.deepcopy(unpad_params)
            dparams[j] += unpad_offsetparams[j]
            J[:,j] = (modelfn(
                get_pad_params(dparams, displ_list, params), merged_list) - base_mod_list)/unpad_offsetparams[j]
    else:
        arg_list = []
        for j in range(len(unpad_params)):
            dparams = copy.deepcopy(unpad_params)
            dparams[j] += unpad_offsetparams[j]
            arg_list.append((get_pad_params(dparams, displ_list, params), merged_list))

        Jtmp = pool.map(modelfn, arg_list)

        Jtmp = [(Jtmp[j] - base_mod_list)/unpad_offsetparams[j] for j in range(len(unpad_params))]
        J = transpose(array(Jtmp))


    if not use_dense_J:
        J = lil_matrix(J)
        J = J.tocsr()

    return J



def get_pad_params(unpad_params, displ_list, params):
    pad_params = copy.deepcopy(params)


    #for i in range(len(displ_list)):
    #    if displ_list[i] == 0:
    #        pad_params = insert(pad_params, i, params[i])
    inds = where(displ_list != 0)
    pad_params[inds] = unpad_params
    return pad_params



def get_unpad_params(pad_params, displ_list):

    unpad_params = copy.deepcopy(pad_params)

    unpad_params = unpad_params.compress(displ_list != 0)

    return unpad_params


def chi2fromresid(resid, Wmat):
    if Wmat == None:
        chi2 = sum(resid**2.)
    else:
        chi2 = dot(dot(resid, Wmat), resid)

    if isnan(chi2):
        print("nan found! ", resid, Wmat)
        return 1.e100
    else:
        return chi2


def miniLM(params, orig_merged_list, displ_list, verbose, maxiter = 150, maxlam = 100000, Wmat = None, jacobian_name = "Jacob.fits", return_wmat = False, use_dense_J = False, pool = None, save_jacobian = True):
    params = array(params, dtype=float64)
    displ_list = array(displ_list, dtype=float64)

    # fix_list -- 1 = fix

    merged_list = copy.deepcopy(orig_merged_list)
    modelfn = orig_merged_list[0]
    del merged_list[0]
    del merged_list[0] #placeholder for inlimit

    lam = 1.e-6
    lamscale = 2.



    converged = 0

    if pool == None:
        curchi2 = modelfn(params, merged_list)
    else:
        curchi2 = modelfn((params, merged_list))

    unpad_offsetparams = get_unpad_params(array(displ_list, dtype=float64), displ_list)*1.e-6


    unpad_params = get_unpad_params(params, displ_list)
    if verbose:
        print("unpad_params ", unpad_params)

    itercount = 0
    was_just_searching = 0
    while lam < maxlam and itercount < maxiter:
        itercount += 1


        if verbose:

            print(len(unpad_offsetparams), len(unpad_params), len(displ_list), len(params), len(curchi2))
        if not was_just_searching:
            Jacob = Jacobian(modelfn, unpad_offsetparams, merged_list, unpad_params, displ_list, params, len(curchi2), use_dense_J, pool = pool)
        #save_img(Jacob.todense(), "tmpjacob.fits")
        was_just_searching = 0

        if verbose:
            print("Jacob.shape ", Jacob.shape)
        Jacobt = transpose(Jacob)

        if verbose:
            print("Dot start ", time.asctime())
        if Wmat == None:
            JtJ = Jacobt.dot(Jacob)
            if not use_dense_J:
                JtJ = JtJ.todense()
        else:
            try:
                JtJ = Jacobt.dot(transpose(Jacobt.dot(Wmat)))
            except:
                print("Couldn't do dot product!")
                print("Sizes ", Jacobt.shape, Wmat.shape)
                sys.exit(1)
        if verbose:
            print("Dot end ", time.asctime())



        JtJ_lam = copy.deepcopy(JtJ)

        for i in range(len(JtJ)):
            JtJ_lam[i,i] *= (1. + lam)
        try:
            if Wmat == None:
                delta1 = -linalg.solve(JtJ_lam, Jacobt.dot(curchi2))
            else:
                delta1 = -linalg.solve(JtJ_lam, Jacobt.dot(dot(Wmat, curchi2)))
        except:

            print("Uninvertible Matrix!")
            if save_jacobian:
                #Jdense = Jacob.todense()
                #print Jdense.shape
                #for i in range(len(Jdense)):
                #    print i
                #    if all(Jdense[i] == 0):
                #        print "fdsakfdskjahfdkjs"
                if not use_dense_J:
                    save_img(Jacob.todense(), jacobian_name)
                else:
                    save_img(Jacob, jacobian_name)

            return [
                array([get_pad_params(unpad_params, displ_list, params)],dtype=float64),
                array([chi2fromresid(curchi2, Wmat), -1], dtype=float64)] + [Jacob.todense()]*return_wmat

        JtJ_lam2 = copy.deepcopy(JtJ)

        for i in range(len(JtJ)):
            JtJ_lam2[i,i] *= (1. + lam/lamscale)
        if Wmat == None:
            delta2 = -linalg.solve(JtJ_lam2, Jacobt.dot(curchi2))
        else:
            delta2 = -linalg.solve(JtJ_lam2, Jacobt.dot(dot(Wmat, curchi2)))

        unpad_params1 = get_pad_params(unpad_params + delta1, displ_list, params)
        unpad_params2 = get_pad_params(unpad_params + delta2, displ_list, params)

        if pool == None:
            chi2_1 = modelfn(unpad_params1, merged_list)
            chi2_2 = modelfn(unpad_params2, merged_list)
        else:
            chi2_1 = modelfn((unpad_params1, merged_list))
            chi2_2 = modelfn((unpad_params2, merged_list))


        if chi2fromresid(chi2_2, Wmat) < chi2fromresid(curchi2, Wmat):
            curchi2 = chi2_2
            unpad_params = unpad_params + delta2
            lam /= lamscale

        elif chi2fromresid(chi2_1, Wmat) < chi2fromresid(curchi2, Wmat):
            curchi2 = chi2_1
            unpad_params = unpad_params + delta1
        else:

            itercount -= 1
            was_just_searching = 1

            while (chi2fromresid(chi2_1, Wmat) >= chi2fromresid(curchi2, Wmat) and lam < maxlam) or isnan(chi2fromresid(chi2_1, Wmat)):

                if verbose:
                    print("Searching... ", lam)
                lam *= lamscale

                JtJ_lam = copy.deepcopy(JtJ)

                for i in range(len(JtJ)):
                    JtJ_lam[i,i] *= (1. + lam)

                if Wmat == None:
                    delta1 = -linalg.solve(JtJ_lam, Jacobt.dot(curchi2))
                else:
                    delta1 = -linalg.solve(JtJ_lam, Jacobt.dot(dot(Wmat, curchi2)))

                unpad_params1 = get_pad_params(unpad_params + delta1, displ_list, params)
                if pool == None:
                    chi2_1 = modelfn(unpad_params1, merged_list)
                else:
                    chi2_1 = modelfn((unpad_params1, merged_list))



        if verbose:
            print("itercount, unpad_params, lam, curchi2 ", itercount, unpad_params, lam, chi2fromresid(curchi2, Wmat))

    if save_jacobian:
        if not use_dense_J:
            save_img(Jacob.todense(), jacobian_name)
        else:
            save_img(Jacob, jacobian_name)

    return [
        array([get_pad_params(unpad_params, displ_list, params)],dtype=float64),
              array([chi2fromresid(curchi2, Wmat)],dtype=float64)] + [JtJ]*return_wmat


def miniLM_new(ministart, miniscale, residfn, passdata, verbose = False, maxiter = 150, maxlam = 100000, Wmat = None, jacobian_name = "Jacob.fits", use_dense_J = False, return_Cmat = True, pad_Cmat = False, pool = None, save_jacobian = False):
    [P, F, param_wmat] = miniLM(ministart, [residfn, None, passdata], miniscale, verbose, maxiter = maxiter, maxlam = maxlam, Wmat = Wmat, jacobian_name = jacobian_name, return_wmat = True, use_dense_J = use_dense_J, pool = pool, save_jacobian = save_jacobian)


    if len(param_wmat) > 0 and return_Cmat:
        if linalg.det(param_wmat) != 0.:
            Cmat = linalg.inv(param_wmat)
        else:
            print("Determinant is zero!")
            Cmat = None
    else:
        Cmat = None
    return P[0], F[0], Cmat

def miniGN(params, orig_merged_list, displ_list, verbose, pool = None):
    params = array(params, dtype=float64)
    displ_list = array(displ_list, dtype=float64)

    # fix_list -- 1 = fix

    merged_list = copy.deepcopy(orig_merged_list)
    modelfn = orig_merged_list[0]
    del merged_list[0]
    del merged_list[0] #placeholder for inlimit

    lam = 1.
    lamscale = 2.


    maxiter = 150

    converged = 0

    curchi2 = modelfn(params, merged_list)

    unpad_offsetparams = get_unpad_params(array(displ_list, dtype=float64), displ_list)*1.e-6


    unpad_params = get_unpad_params(params, displ_list)
    print("unpad_params ", unpad_params)

    itercount = 0
    while itercount < maxiter:
        itercount += 1


        if verbose:
            print(len(unpad_offsetparams), len(unpad_params), len(displ_list), len(params), len(curchi2))

        Jacob = Jacobian(modelfn, unpad_offsetparams, merged_list, unpad_params, displ_list, params, len(curchi2), use_dense_J, pool = pool)

        if verbose:
            print("Jacob.shape ", Jacob.shape)
        Jacobt = transpose(Jacob)
        print("Dot Start")
        JtJ = dot(Jacobt, Jacob)
        print("Dot Finish")


        chi2_1 = 1.e10
        lam = 1.e-10

        JtJ_lam = copy.deepcopy(JtJ)


        try:
            print("delta1")
            delta1 = -dot(linalg.inv(JtJ_lam), dot(Jacobt, curchi2))
        except:

            print("Uninvertible Matrix!")
            return [
                array([get_pad_params(unpad_params, displ_list, params)],dtype=float64),
                array([chi2fromresid(curchi2, Wmat), -1], dtype=float64)]


        while  sum(chi2_1**2.) >= sum(curchi2**2.):


            unpad_params1 = get_pad_params(unpad_params + delta1*lam, displ_list, params)

            chi2_1 = modelfn(unpad_params1, merged_list)

            if sum(chi2_1**2.) < sum(curchi2**2.):
                curchi2 = chi2_1
                unpad_params = unpad_params + delta1*lam
            else:
                lam /= lamscale
                if verbose:
                    print("lam ", lam)



        if verbose:
            print("itercount, unpad_params, lam, curchi2 ", itercount, unpad_params, lam, sum(curchi2**2.))

    return [
        array([get_pad_params(unpad_params, displ_list, params)],dtype=float64),
              array([sum(curchi2**2.)],dtype=float64)]
