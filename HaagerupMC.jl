using Parameters
using PyPlot
using InteractiveUtils


#using Memoize # I was trying to memoize the weight calculations but so far onlw saw the slowdown. Probably the dictionary look-up is as slow as/ slower than calculate it every time. The comment-outs starting from @memoize is the reminicent of this attempt.


@with_kw struct Params
    nu :: Int # Z_nu Haagerup-Izumi category
    zeta ::Float64 = (nu+sqrt(nu^2+4))/2 # qdim(v) for non-invertible v
    N :: Int = 2*nu # number of simples
    A :: Float64 = 1.0 #The unisotropy parameter in AFM model
    height :: Int # Highet of the lattice
    width :: Int # width of the lattice
    niter :: Int
    ival :: Int = nu #initial value
    rec_int :: Int #interval of recordings of the observables
    mdist :: Int = div(height,2) #to what distance the 2pt funcitons are measured
    mpts :: Int = div(width,4) #For each sweep how many points are measured
    n_therm :: Int = 1000
end

function simulate(params, name_dir)
    # variables: 0,1,.. nu-1: invertible ones. nu,nu+1, .. N: non-invertible ones.  N = 2 nu.
    @unpack nu, zeta, N, A, height, width, niter, ival, rec_int, mdist, mpts, n_therm = params
    is_inv(s) = (s <= nu-1) # whther invertible or not
    m :: Int, n :: Int = width, height # just for brevity
    dual(s) = ifelse(is_inv(s),mod(-s,nu),s) # dual of s
    qdim(s) = ifelse(is_inv(s),1,zeta); # quantum dimension
    function W_plaq(x :: Int,z :: Int,y :: Int,yp :: Int)::Float64
        ifelse(x==z,1.0/qdim(x),0.0) + ifelse(y==yp,A/qdim(y),0.0)
    end

    W_ratio_y(x,z,y,yp,new) = W_plaq(x,z,new,yp)/W_plaq(x,z,y,yp) #Ratio when y -> new.
    W_ratio_x(x,z,y,yp,new) = W_plaq(new,z,y,yp)/W_plaq(x,z,y,yp) #Ratio when x -> new.

    #memoize Dict{NTuple{7,Int},Float64} function W_ratio_left(aa::Int ,bb::Int ,ll::Int ,al::Int,bl::Int,p::Int,new::Int)::Float64 
    function W_ratio_left(aa::Int ,bb::Int ,ll::Int ,al::Int,bl::Int,p::Int,new::Int)::Float64 
        W_ratio_x(p,al,aa,ll,new) * W_ratio_y(ll,bb,p,bl,new)
    end

    #@memoize Dict{NTuple{7,Int},Float64} function W_ratio_right(aa::Int ,bb::Int ,rr::Int ,ar::Int,br::Int,p::Int,new::Int)::Float64 
    function W_ratio_right(aa::Int ,bb::Int ,rr::Int ,ar::Int,br::Int,p::Int,new::Int)::Float64 
        W_ratio_x(p,br,bb,rr,new) * W_ratio_y(aa,rr,p,ar,new)
    end

    #@memoize Dict{NTuple{10,Int},Float64} function W_ratio_around_site(aa::Int ,bb::Int ,ll::Int ,rr::Int, al::Int, ar::Int, bl::Int, br::Int, p::Int,new::Int)::Float64
    function W_ratio_around_site(aa::Int ,bb::Int ,ll::Int ,rr::Int, al::Int, ar::Int, bl::Int, br::Int, p::Int,new::Int)::Float64
    #function W_ratio_around_site(aa::Int ,bb::Int ,ll::Int ,rr::Int, al::Int, ar::Int, bl::Int, br::Int, p::Int,new::Int)::Float64
        W_ratio_left(aa,bb,ll,al,bl,p,new)*
        W_ratio_right(aa,bb,rr,ar,br,p,new)*
        qdim(new)/qdim(p)
    end 

    val(ss,(i,j)) = ss[i,j] :: Int
    aa((i,j)) = (mod1(i-1,m),j) #above p
    bb((i,j)) = (mod1(i+1,m),j) #below p
    ll((i,j)) = (i,mod1(j-1,n)) #left p
    rr((i,j)) = (i,mod1(j+1,n)) #right p
    al((i,j)) = (mod1(i-1,m),mod1(j-1,n)) #above left of p
    ar((i,j)) = (mod1(i-1,m),mod1(j+1,n)) #above right of p
    bl((i,j)) = (mod1(i+1,m),mod1(j-1,n)) #below left of p
    br((i,j)) = (mod1(i+1,m),mod1(j+1,n)) #below right of p

    function candidate(v_aa,v_bb,v_ll,v_rr,v_p) #return a candidate form allowed variables (Valid only for the identity projector)
        if v_aa != v_ll || v_aa != v_rr || v_bb != v_ll || v_bb != v_rr
            return v_p
        end
        if v_aa >= nu
            if v_p >= nu
                r = rand(0:nu-1)
                return ifelse(r == 0, mod(-v_aa,nu),ifelse(r+nu-1 < v_p, r+nu-1, r+nu))
            else
                r = rand(nu:2*nu-1)
                return r
            end
        else
            return v_p
        end
    end

    function sweep_once(ss, m :: Int, n :: Int)
        for j in 1:n, i in 1:m
            p = (i,j)
            v_aa,v_p,v_bb,v_ll,v_rr = 
            val(ss,aa(p)),val(ss,p),val(ss,bb(p)),val(ss,ll(p)),val(ss,rr(p))
            cand = candidate(v_aa,v_bb,v_ll,v_rr,v_p)
            if cand != v_p
                v_al,v_bl,v_ar,v_br = 
                val(ss,al(p)),val(ss,bl(p)),val(ss,ar(p)),val(ss,br(p))
                if rand() < W_ratio_around_site(v_aa,v_bb,v_ll,v_rr,v_al,v_ar,v_bl,v_br,v_p,cand)
                    ss[i,j]  = cand
                end
            end
        end
    end
    
    function thermalize_once(ss)
        for j in 1:n, i in 1:m
            p = (i,j)
            v_aa,v_p,v_bb,v_ll,v_rr = 
            val(ss,aa(p)),val(ss,p),val(ss,bb(p)),val(ss,ll(p)),val(ss,rr(p))
            cand = candidate(v_aa,v_bb,v_ll,v_rr,v_p)
            if rand()<(1.0-1.0/N)
                ss[i,j]  = cand
            end
        end
    end

    onepts = [0.0] ::Array{Float64}
    E_2pt = [0.0 for _ in 1:mdist] :: Array{Float64}
    records = [onepts,E_2pt]

    function init_recs!(records)
        fill!(records[1],0.0)
        fill!(records[2],0.0)
    end

    function measure!(ss,records::Array{Array{Float64,1},1}, m :: Int, n :: Int)
        local onepts :: Array{Float64}, E_2pt :: Array{Float64}  = records
        for _ in 1:mpts
            i = rand(1:m) :: Int
            j = rand(1:n) :: Int
            p = (i,j) :: Tuple{Int,Int}
            v_p,v_bb,v_br,v_rr =
            val(ss,p),val(ss,bb(p)),val(ss,br(p)),val(ss,rr(p))
            E_tmp =  - log(W_plaq(v_p,v_br,v_bb,v_rr)sqrt(qdim(v_bb)qdim(v_rr)))
            onepts[1] += E_tmp
            for k in 1:mdist
                pk = (mod1(i+k,m),j) :: Tuple{Int,Int}
                vk_p,vk_bb,vk_br,vk_rr =
                val(ss,pk),val(ss,bb(pk)),val(ss,br(pk)),val(ss,rr(pk))
                E_tmp_k = - log(W_plaq(vk_p,vk_br,vk_bb,vk_rr)sqrt(qdim(vk_bb)qdim(vk_rr)))
                onepts[1] += E_tmp_k 
                E_2pt[k] += E_tmp * E_tmp_k
            end
        end
    end

    mkpath(name_dir)
    onept_name(o)::String = string(name_dir,"/$(o)_1pt.txt") #o is the name of an obervable.
    twopt_name(o,dist)::String = string(name_dir,"/$(o)_$(dist)_2pt.txt") #dist is the distance between the 2pts.


    function write_in_file(records)
        onepts, E_2pt  = records
        open(onept_name("E"), "a+") do io
            print(io,onepts[1]/((1+mdist)*mpts*rec_int)," ")
        end
        for k in 1:mdist
            open(twopt_name("E",k), "a+") do io
                print(io,E_2pt[k]/(mpts*rec_int)," ")
            end
        end
    end

    function thermalize(ss,niter)
        println("thermalize")
        for nn in 1:niter
            thermalize_once(ss)
        end
    end

    function sweep(ss,niter,records,m,n)
        println("start")
        for nn in 1:niter
            sweep_once(ss,m,n)
            if nn % rec_int == 0
                println(nn)
                write_in_file(records)
                init_recs!(records)
            end
            measure!(ss,records,m,n)
        end
    end

    s = [ival for i in 1:m, j in 1:n]
    s0 = deepcopy(s)
    thermalize(s,n_therm)
    s1 = deepcopy(s)
    #@code_warntype W_ratio_around_site(2,2,2,2,2,2,2,2,2,2)
    @time sweep(s,niter,records,m,n)
    #@code_warntype measure!(s,records,m,n)

    #plot
    function plot2d(fig, s)
        m, n = size(s)
        pcolormesh(0:m, 0:n, s, vmin=0, vmax=N-1, cmap="gist_earth")
        fig[:set_aspect]("equal")
    end
    figure(figsize=(12,3)) 
    fig = subplot(141); plot2d(fig, s0); title("Before thermalize")
    fig = subplot(142); plot2d(fig, s1); title("t=0")
    fig = subplot(143); plot2d(fig, s); title("t=$niter")
    tight_layout() 
    savefig(string(name_dir,"/plot.pdf"))
end


function main(args)
    params = Params(nu = parse(Int,args[1]), height = parse(Int,args[2]), width = parse(Int, args[3]), niter = parse(Int, args[4]),rec_int = parse(Int, args[5]))
    name_dir = args[6]
    simulate(params,name_dir)
end

main(ARGS)