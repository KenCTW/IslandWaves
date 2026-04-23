"""
	jLWaves

This set of functions computes island-trapped wave properties in a stratified ocean.
This collection seeks a complex frequency.

To use this collection:

Make an input tuple:

	arr = jLWavesSetup()

Modify an input tuple:

	narr = jLWavesFinch(arr)

To run:

	jLWavesM(arr)


Translated from Matlab bigi___ code  (6/8/2018 version)

K.H. Brink  4/20/2026
"""



"""
	Function jLWavesM(arr)

	Main function for computing stratified island-trapped wave modes
		with complex frequency.

	call as
		jLWavesM(arr)

	When it is done, it will ask if you want to save the results

K.H. Brink 4/20/2026

"""



function jLWavesM(arr)

#   Julia version: file to compute island-trapped wave modes 
#       with stratification, bottom topography, and mean flow.
#       Allows for complex frequency.
#       Modes are not orthogonal.       
#       Calls various other functions beginning with "jLWaves" or other
#
#   Call as 
#       jLWavesM(arr)
#   where
#       arr is an input tuple containing all the pertinent data.
#
#   See function jLWavesSetup to help set this up, and jLWavesFinch to change the input
#       tuple conveniently
#
#   K. Brink   4/20/2026


# arr is a tuple consisting of:
#
#	nn	    number of r grid points
#	mm	    number of vertical grid points
#       wz          where wz is the first guess frequency (rad/sec), entered as an array
#       del         del = 0 for rigid lid; del = 1 for free surface
#       iobc        iobc = 0 says no flow through offshore boundary
#                   iobc = 1 says open offshore boundary
#       f           Coriolis parameter (rad/sec)
#       rmax        offshore size of grid  (km)
#       eps         nominal fractional accuracy for solutions (0.001 works)
#       npts        number of points on dispersion curve
#       nwz         first azimuthal wavenumber (rad/cm)
#       dnw         increment of alongshore wavenumber (rad/cm)
#       ndep        number of offshore distance, depth pairs to read  (in km, m)
#                       Must be >= 1
#       rdep   (ndep times)    radial locations (km) for topography. First value sets the coastline.
#       depr   (ndep times)    depth (m) corresponding to r
#       nR          Number of distance, bottom friction pairs (>=0)
#                       nR = 0 means no bottom friction
#       rR     (nR times)  radial locations for R in km  
#       Rread  (nR times)  bottom resistance coefficient values  in cm/sec
#       nnsq        number of Nsquared values to read  (>=1)
#       zr          depth increment for Nsquared to read  (m). First depth is
#                       at the surface.
#       alph        exponential tail scale (km) for extrapolating Nsquared beyond
#                       the max depth read
#       nsqr    (nnsq times)    Nsquared values (rad^2/s^2)
#       vzero       amplitude of mean azimuthal flow (cm/sec)
#                      If vzero = 0, mean flow is zero everywhere
#       rzero       radial location of maximum in azimuthal mean flow (km)
#       zzero       depth of mean flow maximimum  (m) (> 0 for subsurface max)
#       zscaled     downward e-folding scale for mean velocity (m)
#       zscaleu     upward e-folding scale for mean velocity (m)
#       xscaleoff   offshore radial e-folding scale for mean velocity (km)
#       xscaleon    onshore radial e-folding scale for alongshore velocity (km)
#       kk          determines where Nsquared is undisturbed                  
#        		kk = 1          undisturbed offshore
#                       kk = 0          undisturbed at coast
#       ipause      pause to see graphics (1) or execute without pauses (0)


#       All internal works are in cgs units, although inputs are in convenient units

#       Numerical grid is arranged so there is a row of grid points outside the physical		
#           domain on all boundaries. Thus, there are mm-2 vertical grid points within,			
#           or on the border of, the physical domain.
#       The grid points are arranged so that it is preferable to have nn > mm



	global nn, mm, dt, dr, f, rl, del, iobc, BB, h, rrzero


	#   Extract array size
	nn = arr[1]
	mm = arr[2]
	println("nn, mm = ", nn, ",  ", mm)

	#   First frequency guess
	wgr = arr[3]
	wg = wgr[1] + im*wgr[2]

	#   Rigid lid/free surface, offshore wall or open
	del = arr[4]
	iobc = arr[6] 
	println(' ')

	if del > 0.5
    		println("Free surface")
	else
    		println("Rigid lid")
	end

	if iobc > 0.5
    		println("Open BC at r = rmax")
	else
	        println("Closed BC at r = rmax")
	end

	#   f
	f = arr[8]
	println(' ')
	@printf("f = %.3e rad/sec \n", f)

	#   Domain radial size 
	rmax = arr[9]

	#   Desired fractional accuracy
	eps = arr[10]

	#   number of dispersion curve points, first wavenumber and wavenumber increment
	npts = arr[11]
	nwz = arr[12]
	dnw = arr[13]
	nwl = nwz*ones(npts) + dnw*(0:(npts-1))
	wnn = 0*ones(npts,2)

	#   Extract depth information (#, x, depth)
	ndep = arr[14]
	rdep = arr[15]
	depr = arr[16]

	#   Extract bottom friction information (number, x, r)
	nR  = arr[17]
    	rR = arr[18]
    	Rread = arr[19]

	#   Extract base-state Nsquared information
	nnsq = arr[20]
	zr = arr[21]
	alph = arr[22]
	nsqr = arr[23]

	#   Pull out information on mean alongshore flow
	vzero = arr[24]
    	rzero = arr[25]
    	zzero   = arr[26]
    	zscaled = arr[27]
    	zscaleu = arr[28]
    	rscaleoff = arr[29]
    	rscaleon = arr[30]
    	kk = arr[31]

	#    Get the command for pausing during execution or not
	ipause = arr[32]

#   Now get started: calculate the mean fields

	rmax = rmax*1.0e05
	dt = 1/(mm-3)
	rrzero = rdep[1]*1.0e5
	dr = (rmax-rrzero)/(nn-3)


	icoastflag = jLWavesDep(rdep,depr)
	jSWavesConChk()
	jLWavesRr(nR,rR,Rread)
	jHWavesNsq(zr,alph,nsqr) 

	vtup  = jLWavesVzCal(vzero,rzero,zzero,zscaleu,zscaled,rscaleoff,rscaleon)
	rtup = jLWavesRhoCal(vzero,vtup, kk,ipause)
	
	igravflag = rtup[2]


#   Check to see if there are flags for gravitational instability or inconsistency
	if igravflag != 0
		println(' ')
    		println("Exit because density is gravitationally unstable")
		return NaN
	end

	if icoastflag > 1.5
		println(' ')
    		println("Exit because of attempt to use open r= rmax boundary condition when the bottom is not flat there")
		return NaN
	end

	sleep(1)

#   Set parameters for search
	maxit = 80
	wnn = NaN*ones(npts,1)*(1 + im)
	
	rgr = 0*ones(nn-2,mm-2)
    	zgr = 0*ones(nn-2,mm-2)
    	for n = 1:nn-2
        	rtemp = rrzero + dr*(n-1)
        	rgr[n,:] = rtemp*ones(1,mm-2)/1e5
        	zgr[n,:] = h[n]*(-1*ones(mm-2) + dt*(0:mm-3))/100
    	end
	p = NaN*ones(nn-2,mm-2)*(1+im)
	BB = NaN*ones(nn*mm)*(1 + im)
	exitflag = NaN


#   Now actually do the calculations

	ncal = 1

	while ncal <= npts
		GLMakie.closeall()
    		rl = nwl[ncal]
    		println(' ')
    		@printf("Wavenumber = %.3i ", rl)
    		println(' ')
		println(' ')
 		println("Iterations:")  

		ggg = [real(wg), imag(wg)]
		steppr = 0.025*abs.(real(wg))							#	search simplex
		steppi = 0.025*abs.(imag(wg))	
		ds = [steppr,  steppi]

		
		opt = Opt(:LN_NELDERMEAD,2)
		NLopt.min_objective!(opt,jLWavesDrvr)
		NLopt.xtol_rel!(opt,eps)
		NLopt.maxeval!(opt,maxit)
		NLopt.initial_step!(opt, ds)


		#ggg = [wg]
		(minf, minx, ret) = NLopt.optimize!(opt,ggg)
		
		neval = NLopt.numevals(opt)

		if ret != :XTOL_REACHED								# ret is a "symbol" variable	
			exitflag = 0
			println(' ')
			println("Used up ", maxit, " iterations")
			www = NaN*(1 + im)
		else
			exitflag = 1
			www = minx[1] + minx[2]*im
		end


    
   # output result
   
   	 	if exitflag == 1
        		println(' ')
        		println("Converged!")
        		println(' ')
        		wnn[ncal,1] = www
     
        		if ncal ==1 && npts != 1
            			c = www/rl
            			wg = c*nwl[2]
        		elseif ncal >=2
            			cg = wnn[ncal,1] - wnn[ncal-1,1]
            			wg = www + cg'
            			if ncal > 2.5
                			dc = wnn[ncal,1] - 2*wnn[ncal-1,1] + wnn[ncal-2,1]
                			wg = wg + dc
            			end
        		end
        		#ncal = ncal +1

			fig3 = Figure()
			ax1 = Axis(fig3[1,1])
			scatter!(ax1,nwl,real(wnn[:,1]), marker = :cross, color = :black, markersize = 20)
			scatter!(ax1,nwl,imag(wnn[:,1]), marker = :cross, color = :red, markersize = 20)
			xxx = [0, (maximum(nwl)+1)]
			yyy = [0, 0]
			lines!(ax1,xxx,yyy, color = :black)
			lines!(ax1,nwl,real(wnn[:,1]), color = :black, linestyle = :dash)
			lines!(ax1,nwl,imag(wnn[:,1]), color = :red, linestyle = :dot)
			ax1.ylabel = "ω [1/s]"
			ax1.xlabel = "Azimuthal Wavenumber"
			ax1.title = "Dispersion Curve (Imaginary in red)"
			xlims!(ax1,0.,(maximum(nwl)+1))
			hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )
			display(fig3)
			sleep(1)

#		Calculate and plot pressure
        
       		
        		wj = jLWavesDrvr(minx,0.)
			csign = conj(BB[(2*mm-1)])/abs2(BB[(2*mm -1)])
        		
    			BB = BB*csign

    			for n = 1:nn-2
        			mlow = n*mm + 2
        			mhigh = (n+1)*mm - 1
        			p[n,:] = BB[mlow:mhigh]
    			end
			p = jLWavesNorm(p)


			if ipause > 0.5				

				rff = rgr[:,1]
				hh = -zgr[:,1]
				ybb = -maximum(h)*ones(nn-2)
    
				fig4 = Figure()
				ax1 = Axis(fig4[1,1])

				tol = 0.001
				rmm = maximum(maximum(abs.(real(p))))
				amm = maximum(maximum(abs.(imag(p))))
				rata = amm/rmm
				pfac = 1
				if rata <= tol
					pfac = NaN
				end

				civ = jHWavesCCon(p,10)

				contour!(ax1, rgr,zgr,real(p), levels =civ, color = :black)	
				contour!(ax1, rgr,zgr,real(p),levels = vec([0 0]), linewidth = 3, color = :black)
				contour!(ax1, rgr,zgr,pfac*imag(p), levels = civ, color = :red)	
				contour!(ax1, rgr,zgr,pfac*imag(p),levels = vec([0 0]), linewidth = 3, color = :red)


				band!(ax1,rff,ybb,-hh,color = :blue)
				xlims!(ax1,minimum(rgr[:,1]),maximum(rgr[:,1]))
				ylims!(ax1,minimum(zgr[:,1]),0)
				
				ax1.xlabel = "x [km]"
				ax1.ylabel = "z [m]"
				ax1.title = "Pressure (heavy contour is 0, real is in black, imaginary is in red)"
				hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )
				
				display(GLMakie.Screen(),fig4)

				if ncal < npts 
					println(" ")
					println("Sleeping for 10 seconds")
					sleep(10)
				end	
			end
			ncal = ncal + 1
			
		else
			ncal = npts + 1
			
		end
		
	end
						# end of ncal loop

	rtup = jLWavesRhoCal(vzero,vtup, kk,0)


	
	if ipause  < 0.5 && exitflag == 1
		p = jLWavesNorm(p)

		rff = rgr[:,1]
		hh = -zgr[:,1]
		ybb = -maximum(h)*ones(nn-2)
    
		fig4 = Figure()
		ax1 = Axis(fig4[1,1])

		tol = 0.001
		rmm = maximum(maximum(abs.(real(p))))
		amm = maximum(maximum(abs.(imag(p))))
		rata = amm/rmm
		pfac = 1
		if rata <= tol
			pfac = NaN
		end

		cip = jHWavesCCon(p,10)		

		contour!(ax1, rgr,zgr,real(p), levels =cip, color = :black)	
		contour!(ax1, rgr,zgr,real(p),levels = vec([0 0]), linewidth = 3, color = :black)
		contour!(ax1, rgr,zgr,pfac*imag(p), levels = cip, color = :red)	
		contour!(ax1, rgr,zgr,pfac*imag(p),levels = vec([0 0]), linewidth = 3, color = :red)

		band!(ax1,rff,ybb,-hh,color = :blue)
		xlims!(ax1,minimum(rgr[:,1]),maximum(rgr[:,1]))
		ylims!(ax1,minimum(zgr[:,1]),0)
				
		ax1.xlabel = "r [km]"
		ax1.ylabel = "z [m]"
		ax1.title = "Pressure (heavy contour is 0, real is in black, imaginary is in red)"
		hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )
				
		display(GLMakie.Screen(),fig4)
	end			

	println(' ')
	if exitflag  == 1
		jLWavesEngDiag(wnn[ncal-1,1],BB)
	end


#	Plot final dispersion curve (real and complex)

	if npts > 1.5 
		fig2 = Figure()
		ax1 = Axis(fig2[1,1])
		#ax2 = Axis(fig2[2,1])
		scatter!(ax1,nwl,real(wnn[:,1]), marker = :cross, color = :black, markersize = 20)
		scatter!(ax1,nwl,imag(wnn[:,1]), marker = :cross, color = :red, markersize = 20)

		lines!(ax1,nwl,real(wnn[:,1]), color = :black, linestyle = :dash)
		lines!(ax1,nwl,imag(wnn[:,1]), color = :red, linestyle = :dot)
		xxx = [0, (maximum(nwl)+1)]
		yyy = [0, 0]
		lines!(ax1,xxx,yyy, color = :black)


		ax1.ylabel = "ω [1/s]"
		ax1.xlabel = "Azimuithal Wavenumber"
		ax1.title = "Dispersion Curve (Real in black, Imaginary in red)"
		xlims!(ax1,0.,(maximum(nwl)+1))
		hidedecorations!(ax1; label = false, ticklabels = false, ticks = false )

		display(fig2)
	
	end

    
#       Option for saving results

	isave = inputi("Do you want to save this information? 0 (no) or 1 (yes) ")
    	if isave > 0.5
		println("Enter the name of the desired file (.jld2 will be appended to this) ")	
		str = readline()
		fname = str*".jld2"
		ilong = NaN
		otup = [f, iobc, del, nwl, wnn, rgr, zgr, p, rR, Rr, vzero]
		save_object(fname,otup)
		
		println(' ')
		println("Saved as a tuple in file ", fname)
		verbiage = """
			Contents of the tuple:

				1) f	(Coriolis parameter: rad/sec)
				2) iobc	 (= 0 for wall at offshore boundary, = 1 for open)
				3) del	 (= 0 for rigid lid, = 1 for free surface)
				4) nwl (wavenumbers)
				5) wnn (complex frequencies in rad/sec)
				6) rgr (r grid in km)
				7) zgr (z grid in m)
				8) p  (pressure)
				9) rR (R for bottom frictional coefficient: km)
				10) Rr (bottom frictional coefficient: cm/sec)
				11) vzero  (extreme value of mean velocity: cm/sec)

			You can load this tuple with the 'load_object' command.
			"""
			
		println(verbiage) 
         
			
    	end

end				#  end of function




"""	
	jLWavesDep(rdep,depr)

	Calculate depth and its derivative for an island case

K.H. Brink 4/1/2026
"""


function jLWavesDep(rdep,depr)

#       Compute the depth profile, given the inputs
#       Only compute depth for points within the domain or on the edges
#       iflag != 0 signals trouble with using open offshore boundary
#           i.e. hr != 0    
#
#       K.H. Brink 4/1/2026

	global nn, dr, h, hr, hrr, iobc, rrzero


	r = rrzero .+ dr*(0:(nn-3))
	rmax = r[nn-2]

	h = 0*ones(nn-2)
	hr = 0*copy(h)
	hrr = 0*copy(h)

	rdep = rdep*1.0e5
	depr = depr*100.

	if rmax > rdep[end]
		rdep[end] = rmax
	end

	n = length(rdep)
    	dxs = dr*dr
    	dx2 = 2*dr

	if n == 1
    		h = depr*ones(nn-2)
    		hr = 0*h
    		hrr = hr
	else


		itp = interpolate((vec(rdep),),vec(depr),Gridded(Linear()))
		h = itp(r)

 		hr[2:(nn-3)] = (h[3:(nn-2)]-h[1:(nn-4)])/dx2
    		hrr[2:(nn-3)] = (h[3:(nn-2)] - 2*h[2:(nn-3)] + h[1:(nn-4)])/dxs


    		hr[1] = (h[2]-depr[1])/dx2
    		hrr[1] = (depr[1] -2*h[1] + h[2])/dxs
    
    		hr[nn-2] = (h[nn-2]-h[nn-3])/dx2
    		hrr[nn-2] = (h[nn-3] - h[nn-2])/dxs
    
	end

	iflag = 0

	if iobc > 0.5
    		if hr[nn-2] !=0
        		iflag = 2
    		end
	end
	return iflag

end						# end of function










"""
	jSWavesConChk


	Check for consistency of bottom slope
		in cylindrical coordinates

K.H. Brink  11/21/2025
"""

function jSWavesConChk()

#       Check to find max of (hr/h)*(dr/dt)
# 		in cylindrical coordinates
#    K.H. Brink 11/21/2025 based on Matlab version of  3/8/2004

	global hr, h, dr, dt

	arrt = (hr./(h + 0.5*hr*dr))*dr/dt
	arrt = abs.(arrt)

	rmax, ii = findmax(vec(arrt))

	println(' ')
	xm = dr*(ii-1)/1e5
	@printf(" Max consistency ratio = %.3f at r = %.2f km \n", rmax, xm)
	println("         This should be kept less than one, and definitely less than 10")
	println(' ')
end										# end of function





"""
	jLWavesRr(nR,rR,Rread)

	Interpolate bottom friction coefficient onto the grid
		for an island


K.H. Brink   4/1/2026
"""

function jLWavesRr(nR,rR,Rread)

#   Compute the array of bottom friction parameters, given the inputs
#		for an island
#
#	K.H. Brink 4/1/2026

	global nn, dr, R, Rr,rrzero

	if nR == 0
      		R = 0*ones(nn-2)
      		Rr = copy(R)
  	elseif nR ==1
      		R = Rread[1]*ones(nn-2)
      		Rr = 0*ones(nn-2)
  	else

    		rad = rrzero .+ dr*(0:(nn-3))
    		rmax = rad[nn-2]
    		rR = rR*1.0e5
		if rmax > rR[end]
			rR[end] = rmax
		end



		R = 0*ones(nn-2)
      		Rr = 0*ones(nn-2)	

    		rRmax = rR[end]

		itp = interpolate((vec(rR),),vec(Rread),Gridded(Linear()))
		R = itp(rad)

    		tdx = 2*dr

    		Rr[2:(nn-3)] = (R[3:nn-2]- R[1:(nn-4)])/tdx

    		Rr[nn-2] = 0
    		Rr[1] = (R[2]-R[1])/dr

	end
end									# end of function







"""
	jHWavesNsq(zr,alph, nsqr)

	Compute undisturbed N^2 profile


K.H. Brink 11/21/2025
"""




function jHWavesNsq(zr,alph,nsqr)

#  Compute the undisturbed (no mean flow) N^2(z) profile
#       given the inputs
#
#   Gives back nsq(mm-2), the array of interior N^2 points, starting at the bottom
#       (consistent with grid point numbering)
#
#   K.H. Brink 11/21/2025 based on Matlab 5/5/2019 version

	global  mm, dt, h, nsq

	zr = zr*100
	alph = alph*1e5
	n = length(nsqr)

	nsq = ones(mm-2)
    	zmax = maximum(h)
    	nnmax = Int(ceil(zmax/zr) +2)
   	nsqtemp = 0*ones(nnmax)
    	zf = -zr*(0:(nnmax-1))
    	zrr = -zr*(0:(n-1))

    	if zmax > maximum(-zrr)
        	nsqtemp[1:n] = nsqr
        	ii = findall(zf -> zf < minimum(zrr),zf)
        	zfill = zf[ii]
        	nsqtemp[(n+1):nnmax] = nsqr[n]*exp.((zfill .- minimum(zrr))/alph)
        	zzz = zf
    	else
        	zzz = zrr
        	nsqtemp = nsqr
    	end
    	zint = -dt*zmax*(0:(mm-3))

	zzzt = reverse(zzz)
	nsqtempt = reverse(nsqtemp)
    
	itp = interpolate((vec(zzzt),),vec(nsqtempt), Gridded(Linear()))
	nsq = itp(zint)
   
    	nsq = nsq[(mm-2):-1:1]

     	ij = findall(nsq -> 0*nsq !=0, nsq)
    	if isempty(ij) == 0
        	zzzz = min(nsq)
        	nsq[ij] = zzzz
    	end
end							# end of function





""" 
	jLWavesVzCal(vzero,rzero,zzero,zscaleu,zscaled,rscaleoff,rscaleon)

	Compute the mean velocity field for an island case

K.H. Brink  11/21/2025
"""


function  jLWavesVzCal(vzero,rzero,zzero,zscaleu,zscaled,rscaleoff,rscaleon)

#   compute mean velocity field, given inputs
#	for an island case
#   vmean,vmr,vmrr are in sigma coordinates
#   vtemp, vztemp in z coordinates
#       vtemp, vztemp to have size(nn-2,mm-2)    
# 
#       K.H. Brink 4/1/2026

	global nn, mm, dt, dr, h, vm, vmr, vmrr, rrzero

	rzero = rzero*1.0e5
	rscaleoff = rscaleoff*1.0e5
	rscaleon = rscaleon*1.0e5
	zzero = -zzero*100
	zscaleu = zscaleu*100
	zscaled = zscaled*100

	vm = 0*ones(nn-2,mm-2)
	vmr = 0*ones(nn-2,mm-2)
	vmrr = 0*ones(nn-2,mm-2)
	vtemp = 0*ones(nn-2,mm-2)
	vztemp = 0*ones(nn-2,mm-2)

	if abs(vzero) > 0.0001
#       Compute vmean on rectanguar grid first (vtemp, vztemp), then
#           interpolate into sigma (vm)


    		rfactoff = rscaleoff^2
    		rfacton = rscaleon^2
    		zfactu = zscaleu^2
    		zfactd = zscaled^2
    
    		hmax = maximum(h)
    		vtemp = NaN*ones(nn-2,mm-2)
    		r = rrzero .+ dr*(0:(nn-3))
    		z = hmax*(-1*ones(mm-2) + dt*(0:(mm-3)))
    		dz = hmax*dt
    
    		ii = findall(r ->  r < rzero, r)
    		iic = findall(r -> r >= rzero, r)
    		jj = findall(z -> z < zzero, z)
    		jjc = findall(z -> z >= zzero, z)
		slop = 1.0/rzero
    
   		if isempty(ii) != true
        		for n = 1:maximum(ii)
            			if isempty(jj) != true
					nzz = length(jj)
					zzn = zzero*ones(nzz)
                			vtemp[n,jj] = exp.(-((r[n]-rzero).*(r[n]-rzero)/rfacton))*exp.(-(z[jj]-zzn).*(z[jj]-zzn)/zfactd)
            			end
            
            			if isempty(jjc) != true
					nzz = length(jjc)
					zzn = zzero*ones(nzz)
                			vtemp[n,jjc] = exp.(-((r[n]-rzero).*(r[n]-rzero)/rfacton))*exp.(-(z[jjc]-zzn).*(z[jjc]-zzn)/zfactu)
            			end
				#vtemp[n,:] = vtemp[n,:]*r[n]*slop
        		end    
    		end
        
    		if isempty(iic) != true
        		for n = minimum(iic):nn-2
            			if isempty(jj) != true
					nzz = length(jj)
					zzn = zzero*ones(nzz)
                 			vtemp[n,jj] = exp.(-((r[n]-rzero).*(r[n]-rzero)/rfactoff))*exp.(-(z[jj]-zzn).*(z[jj]-zzn)/zfactd)
            			end
            
            			if isempty(jjc) != true
					nzz = length(jjc)
					zzn = zzero*ones(nzz)
                	 		vtemp[n,jjc] = exp.(-((r[n]-rzero).*(r[n]-rzero)/rfactoff))*exp.(-(z[jjc]-zzn).*(z[jjc]-zzn)/zfactu)
            			end          
        		end
        
    
    		end
		vtemp = vzero*vtemp
   
    		dz2 = 2*dz    	
		dr2 = 2*dr
    		drs = dr*dr

    		vztemp[:,2:(mm-3)] = (vtemp[:,3:mm-2]-vtemp[:,1:(mm-4)])/dz2
    		vztemp[:,1] = (vtemp[:,2]-vtemp[:,1])/dz
    		vztemp[:,(mm-2)] = (vtemp[:,(mm-2)]-vtemp[:,(mm-3)])/dz
  
    		vmc = vtemp
    		vmrc = NaN*ones(nn-2,mm-2)
    		vmrrc = NaN*ones(nn-2,mm-2)
    
    		for n = 2:nn-3
        		vmrc[n,:] = (vmc[(n+1),:] - vmc[(n-1),:])/dr2
        		vmrrc[n,:] = (vmc[(n+1),:] - 2*vmc[n,:] + vmc[(n-1),:])/drs
   		end
    
    		vmrc[1,:] = (vmc[2,:]-vmc[1,:])/dr
    		vmrrc[1,:] = 0*vmc[1,1:mm-2]
    
    		vmrc[nn-2,:] = (vmc[nn-2,:] - vmc[nn-3,:])/dr
    		vmrrc[nn-2,:] = 0*vmc[nn-2,1:mm-2]
   	
       

    		vm = jHWavesZtoSig(vmc)				# in sigma coordinates
    		vmr = jHWavesZtoSig(vmrc)
    		vmrr = jHWavesZtoSig(vmrrc)  

	end 

	vtup = (vtemp, vztemp)
	return vtup

end						# end of function







"""
	jLWavesRhoCal(vzero, vtup, kk, ipause)

	Calculate the density field etc., accounting for the mean flow
		for an island case


K.H. Brink 4/1/2026
"""

function jLWavesRhoCal(vzero,vtup,kk,ipause)

#   Compute the overall density field, given the basic N^2 and mean v
#		for an island
#   Then compute N^2, M^2 etc.
#	call as
#		rtup = jLWavesRhoCal(vzero, vtup, kk, ipause)
#
#   Does calculations in z coordinates and then converts to sigma
#
#   Arrays such as n2 have dimension (nn-2,mm-2)       
#  
#	K.H. Brink,  4/1/2026

	global nn, mm, dt, dr, f, h, hr, n2, n2z, m2, m2z, m2r, nsq, vm, vmr, rrzero

	rho = 0*ones((nn-2),(mm-2))

	vtemp = vtup[1]
	vztemp = vtup[2]
	r = rrzero .+ dr*(0:(nn-3))
	
	conts = -1.03/980
	f2 = copy(0*rho)
	f3 = 0*ones((nn-2),(mm-2))
	for ii = 1:nn-2
		f2[ii,:] = f .+ vmr[ii,:] + vm[ii,:]/r[ii]		
    		f3[ii,:] = f .+ 2.0*vtemp[ii,:]/r[ii]			
	end
	cc = -(f3*1.03/980)
	
	hmax = maximum(h)
	rmax = r[end]/1e5		
	rhsave = 0.
	dz = hmax*dt

#   Compute background density  (as if vm = 0)(in Cartsian coordinates)
	rhob = 0*ones(mm-2)
	rhob[1:mm-2] = dz*conts*(cumsum(nsq[1:(mm-2)]) - 0.5*(nsq[1]*ones(mm-2) + nsq[1:(mm-2)]))
	rhob = rhob - minimum(rhob)*ones(mm-2)
	rhop = 0*ones(nn-2,mm-2)

	for n = 1:nn-2
    		rho[n,:] = rhob
	end

	if abs(vzero) > 0.001
    		for m = 1:(mm-2)
        		rhop[1,m] = 0.       
        		#rhop[2:(nn-2),m] = dr*cc[2:nn-2,m].*(cumsum(vztemp[2:(nn-2),m]) + 0.5*( vztemp[1,m]*ones(nn-3) - vztemp[2:nn-2,m]))
        		rhop[2:(nn-2),m] = dr*(cumsum(cc[2:(nn-2),m].*vztemp[2:(nn-2),m]) + 0.5*( cc[1,m]*vztemp[1,m]*ones(nn-3) - cc[2:nn-2,m].*vztemp[2:nn-2,m]))

			if m == 1
            			(yyy,ii) = findmax(vec(abs.(rhop[:,m])))
            			rhsave = rhop[ii,m]
        		end
        		rhop[:,m] = rhop[:,m] - kk*rhop[nn-2,m]*ones(nn-2)
        		rho[:,m] = rho[:,m] + rhop[:,m]
    		end


    #   This assures that the maximum density change is not at the very bottomi
    		rho = rho - rhsave*ones((nn-2),(mm-2))
 
	
	end



	ccc = -980/1.03
	dz2 = 2*dz
	dzs = dz*dz
	dr2 = 2*dr
	drs = dr*dr
	drdt = dr*dz
	drdt4 = drdt*4

	n2c = 0*ones(nn-2,mm-2)
	n2zc = 0*ones(nn-2,mm-2)
	m2c = 0*ones(nn-2,mm-2)
	m2rc = 0*ones(nn-2,mm-2)
	m2zc = 0*ones(nn-2,mm-2)

	for m = 2:mm-3
       		n2c[:,m] = (rho[:,m+1]-rho[:, m-1])/dz2
       		n2zc[:,m] = (rho[:,m+1] - 2*rho[:,m] + rho[:,m-1])/dzs
	end

    	n2c[:,mm-2] = (rho[:,mm-2]-rho[:,mm-3])/dz
    	n2zc[:,mm-2] = 0*rho[:,mm-2]
    
    	n2c[:,1] = (rho[:,2]-rho[:,1])/dz
    	n2zc[:,1] = 0*rho[:,1]
    
    	n2c = n2c*ccc
    	n2zc = n2zc*ccc
    
	for n = 2:nn-3
    		m2c[n,:] = (rho[n+1,:]-rho[n-1,:])/dr2
    		m2rc[n,:] = (rho[n+1,:] - 2*rho[n,:] + rho[n-1,:])/drs
    		for m = 2:mm-3
        		m2zc[n,m] = (rho[n+1,m+1] - rho[n+1,m-1] - rho[n-1,m+1] + rho[n-1,m-1])/drdt4
    		end
    		m2zc[n,1] = 2*(rho[n+1,2] - rho[n+1,1] - rho[n-1,2] + rho[n-1,1])/drdt4
    		m2zc[n,mm-2] = 2*(rho[n+1,mm-2] - rho[n+1,mm-3] - rho[n-1,mm-2] + rho[n-1,mm-3])/drdt4
	end
  
    	m2c[1,:] = (rho[2,:] - rho[1,:])/dr
    	m2rc[1,:] = 0*copy(rho[1,:])
    	m2zc[1,2:mm-3] = 2*(rho[2,3:mm-2] - rho[2,1:mm-4] - rho[1,3:mm-2] + rho[1,1:mm-4])/drdt4
    	m2zc[1,1] = (rho[2,2] - rho[1,2] - rho[2,1] + rho[1,1])/drdt
    	m2zc[1,mm-2] = (rho[2,mm-2] - rho[2,mm-3] - rho[1,mm-2] + rho[1,mm-3])/drdt
	
   	m2c[nn-2,:] = (rho[nn-2,:] - rho[nn-3,:])/dr
    	m2rc[nn-2,:] = 0*copy(rho[nn-2,:])
    	m2zc[nn-2,2:mm-3] = 2*(rho[nn-2,3:mm-2] - rho[nn-2,1:mm-4] - rho[nn-3,3:mm-2] + rho[nn-3,1:mm-4])/drdt4
    	m2zc[nn-2,1] = (rho[nn-2,2] - rho[nn-3,2] - rho[nn-2,1] + rho[nn-3,1])/drdt
    	m2zc[nn-2,mm-2] = (rho[nn-2,mm-2] - rho[nn-2,mm-3] - rho[nn-3,mm-2] + rho[nn-3,mm-3])/drdt

    	m2c = ccc*m2c
    	m2rc = ccc*m2rc
    	m2zc = ccc*m2zc
    
# convert to sigma coordinates

	m2 = jHWavesZtoSig(m2c)
	m2r = jHWavesZtoSig(m2rc)
	m2z = jHWavesZtoSig(m2zc)
	n2 = jHWavesZtoSig(n2c)
	n2z = jHWavesZtoSig(n2zc)


	fs = f2 - (m2.*m2)./(f*n2)
	fstest = minimum(minimum(fs/f))

	ii = findall(n2 -> n2 < 0, n2)


#  	     Make plots of density, mean v

	rgr = 0*ones(nn-2,mm-2)
	zgr = 0*copy(rgr)
	zzgr = 0*copy(zgr)
	for n = 1:nn-2
    		rtemp = (rrzero + dr*(n-1))/1e5
    		rgr[n,:] = rtemp*ones(1,mm-2)
    		zgr[n,:] = h[n]*(-1*ones(mm-2) + dt*(0:mm-3))/100
		zzgr[n,:] = hmax*(-1*ones(mm-2) + dt*(0:mm-3))/100
	end

	zd = maximum(h)
	z = zd*(-1*ones(mm-2) + dt*(0:mm-3))/100
	r = (rrzero .+ dr*(0:nn-3))*1e-5
	hh = h/100

	nqq = length(h)
	rff = (rrzero .+ dr*(0:(nqq-1)))*1e-5
	ybb = -maximum(hh)*ones(nqq)
	rffm = maximum(rff)

	fig2 = Figure()
	ax1 = Axis(fig2[1,1])	
	ax2 = Axis(fig2[2,1])
	

	contour!(ax1, r,z,1000*rho, levels = 10, color = :black)
	
	if isempty( ii) == false
		scatter!(ax1,vec(rgr[ii]),vec(zgr[ii]), color = :red, marker = :cross)
	end
	if (fstest/f) < 0
	#		fs changes sign somewhere!
		contour!(ax1,rgr,zgr,fs*1e5,color = :red, levels = 5)
		println(' ')
		println("fs/f becomes negative somewhere: symmetric instability is possible")
		println("	f* X 10^5 is contoured in red along with density")
		println(' ')
		ax1.title = "Density (sigma-t units) and f* X10^5 (1/sec)"
	else
		ax1.title = "Density (sigma-t units)"
	end
	band!(ax1,rff,ybb,-hh,color = :blue)
	ax1.xlabel = "r [km]"
	ax1.ylabel = "z [m]"
	xlims!(ax1,rrzero/1e5, rmax)
	ylims!(ax1,-zd/100, 0)
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)




	band!(ax2,rff,ybb,-hh,color = :blue)
	if abs(vzero) < 0.001
		text!(ax2,0.6*rffm,-zd*0.6/100, text = "v_0 = 0")
	else
		contour!(ax2,rgr,zzgr,vtemp, color = :black, levels = 5, labels = true, labelsize = 12)
	end
	ax2.xlabel = "r [km]"
	ax2.ylabel = "z [m]"
	tstr = @sprintf("Mean Alongshore velocity (vzero = %.1f cm/sec)", vzero)
	ax2.title = tstr
	xlims!(ax2,rrzero/1e5, rmax)
	ylims!(ax2,-zd/100, 0)

	band!(ax2,rff,ybb,-hh,color = :blue)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)

	
	display(GLMakie.Screen(),fig2)


	if ipause > 0.5
		println(' ')
		println("Sleeping for 10 seconds")
		println(' ')
		sleep(10)
	else
		sleep(0.1)
	end
	

	iflag = 0
	if isempty(ii) == false
    		iflag = 1 
	end

	rtup = (rho, iflag)
	return rtup
end								# end of function






"""
	jHWavesZtoSig

	Convert a field from z to sigma coordinates


K. Brink 11/25/2025
"""



function jHWavesZtoSig(gridz)

#   convert gridz (on x,z grid) to gridsig (on x,sigma grid)
#   works for arrays like vm that are (nn-2,mm-2) 
#   Call as
#		gridisg = jHWavesZtoSig(gridz)
#
#  K.H. Brink 11/25/2025  based on Matlab of 5/24/2019

	global nn, mm, dt, h

	th = -1*ones(mm-2)  + dt*(0:(mm-3))	
	zin = maximum(h)*th
	gridsig = NaN*ones(nn-2,mm-2)

	for n = 1:(nn-2)
    		z = th*h[n]
    		zzz = zin
    		if z[end] > zzz[end]
       		 	zzz[end] = z[end]
    		end
		itp = interpolate((vec(zzz),),vec(gridz[n,:]), Gridded(Linear()))
		gridsig[n,:] = itp(z)
	end
	return gridsig
end							# end of function






"""
	jHWavewsSigtoZ(gridsig)

	Convert sigma to z coordinates

K.H. Brink 12/9/2025
"""



function jHWavesSigtoZ(gridsig)

#   Convert gridded field to z coordinates from sigma coordinates
#	call as
#		gridz = jWavesSigtoZ(gridsig)
#
#	Not used, but provided for completeness
#	K.H. Brink 12/9/2025, based on Matlab of 5/5/03

	global nn, mm, dt, h

	th = -1*ones(mm-2)  + dt*(0:mm-3)
	zout = maximum(h)*th

	gridz = NaN*ones(nn-2,mm-2)

	for n = 1:nn-2
    		z = h[n]*th
    		jj = findall(zout -> zout >= -h(n) && zout <= 0, zout)
		itp = interpolate((vec(z),),vec(gridsig[n,:]),Gridded(Linear()))
		gridz[n,jj] = itp(zout[jj])

    		nss = maximum(size(jj))
    		if nss <(mm-2)
        		jj = find((-1.2*h[n]) < zout &  zout < -h[n]|(zout < -h[n] && zout > -0.2*maximum(h)))
        		gridz[n,jj] = gridsig[n,1]*ones(size(gridz[n,jj]))
    		end  
	end
	return gridz
end								# end of function



"""
		jLWavesDrvr(ww,yyyy)

	set up and solve the matrix equation for an island

K. Brink  4/1/2026
"""

function jLWavesDrvr(ww,yyyy)

#   Create the big array that gets solved for pressure
#		for an island
#   BB is the pressure field in response to the arbitrary forcing at complex frequency w
#		yyyy is a dummy
#   K. Brink 4/1/2026, based on Matlab of 6/27/2018

	global nn, mm, f, dr, h, BB, n2

	w = ww[1] + im*ww[2]
	fsq = f*f
	wf2 = w*fsq
	nm = nn*mm
	nrank = nn*mm



#   Do the corners: multiple choices are possible, but these seem to work

#       1,1 (lower, onshore)  pz = 0

    	aa1 = 1.
 	aa2 = -2.
	aa3 = 1.
	aa4 = [aa1,aa2,aa3]*wf2/n2[1,1]
	AA = sparse([1,1,1],[1,2,3],aa4,nm,nm)


#  	upper, onshore  pz = 0

	aa1 = 1.
    	aa2 = 0
	aa3 = -1.
	aa4 = [aa3, aa2, aa1]*wf2/n2[1,mm-2]
	AA = AA + sparse([mm,mm,mm],[(mm-2), (mm-1),(mm)],aa4,nm,nm)

    
#       Offshore, bottom   pzz = 0

    	mr = (nn-1)*mm +1
	aa1 = 1.
	aa2 =-2.
	aa3 = 1.
	aa4 = [aa1, aa2, aa3]*wf2/n2[nn-2,1]
	AA = AA + sparse([mr, mr, mr],[mr, mr+1, mr+2],aa4,nm,nm)


#       Offshore, top  pz = 0

	aa1 = 1.
	aa2 = 0.
	aa3 = -1.
	aa4 = [aa3, aa2, aa1]*wf2/n2[nn-2,mm-2]
	AA = AA + sparse([nm, nm, nm],[nm-2, nm-1, nm],aa4,nm,nm)

    
    
#   Do the coastal boundary

	cc = jLWavesCc(w)

	jl = 2
	jh = mm-1
	ipts = jl:jh

	vv = vec(cc[:,2])
	AA = AA + sparse(vec(ipts),vec(ipts),vv,nm,nm)

	id = ipts -1*ones(mm-2)						# column
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,3]),nm,nm)

	id = ipts + (mm -1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,4]),nm,nm)

	id = ipts + mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,5]),nm,nm)

	id = ipts + (mm + 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,6]),nm,nm)

	id = ipts + (2*mm - 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,7]),nm,nm)

	id = ipts + 2*mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,8]),nm,nm)

	id = ipts + (2*mm + 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)


  	ibot = (mm+2)*ones(mm-2)  
	AA = AA + sparse(vec(ipts),vec(ibot),vec(cc[:,12]),nm,nm)

	id = ibot - 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,11]),nm,nm)

	id = ibot + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,13]),nm,nm)

	id = ibot -mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,10]),nm,nm)

	id = ibot + mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,14]),nm,nm)

	cc = nothing



#   Do the outer boundary
#   This assumes the bottom is flat here when BC is open

	cc = jLWavesOff(w)
	nlow = nm - mm + 2
	nhigh = nm -1
	ipts =  nlow:nhigh

	AA = AA + sparse(vec(ipts),vec(ipts),vec(cc[:,8]),nm,nm)

	id = ipts -(2*mm +1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts - 2*mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

	id = ipts -(2*mm - 1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,3]),nm,nm)

	id = ipts -(mm +1)*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,4]),nm,nm)

	id = ipts -mm*ones(mm-2)
	AA = AA + sparse(vec(ipts), vec(id), vec(cc[:,5]),nm,nm)

	id = ipts -(mm - 1)*ones(mm-2)	
	AA = AA + sparse(vec(ipts), vec(id), vec(cc[:,6]),nm,nm)

	id = ipts - 1*ones(mm-2)
	AA = AA + sparse(vec(ipts), vec(id), vec(cc[:,7]),nm,nm)

	id = ipts + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)


	ibot = (nm-2*mm + 2)*ones(mm-2)

	AA = AA + sparse(vec(ipts),vec(ibot),vec(cc[:,12]),nm,nm)
	
	id = ibot -1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,11]),nm,nm)
	
	id = ibot + 1*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,13]),nm,nm)
	
	id = ibot -mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,10]),nm,nm)

	id = ibot + mm*ones(mm-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,14]),nm,nm)
	
	cc = nothing





#   Do the bottom boundary condition


	cc = jLWavesBot(w)

	ipts = (mm+1):mm:(nrank-2*mm+1)

	AA = AA + sparse(vec(ipts),vec(ipts),vec(cc[:,4]),nm,nm)

	id = ipts - mm*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts -(mm - 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

	id = ipts -(mm - 2)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,3]),nm,nm)

	id = ipts + 1*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,5]),nm,nm)

	id = ipts + 2*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,6]),nm,nm)

	id = ipts + mm*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,7]),nm,nm)

	id = ipts + (mm + 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,8]),nm,nm)

	id = ipts + (mm + 2)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)


	cc = nothing



#   Do the surface boundary condition


	cc = jLWavesSbc(w)

	ipts = (2*mm):mm:(nrank-mm)

	AA = AA + sparse(vec(ipts),vec(ipts) , vec(cc[:,4]),nm,nm)

	id = ipts - (mm + 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,1]),nm,nm)

	id = ipts - 2*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

	id = ipts - 1*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,3]),nm,nm)

	id = ipts + (mm - 1)*ones(nn-2)
	AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,5]),nm,nm)

	cc = nothing




#   Do the interior points

	for n = 2:(nn-1)

    		cc = jLWavesMatco(w,n)
    		ipts = ((n-1)*mm + 2):1:((n*mm)-1)


		AA = AA + sparse(vec(ipts),vec(ipts),vec(cc[:,5]),nm,nm)

		id = ipts - (mm+1)*ones(mm-2)
		AA = AA + sparse(vec(ipts), vec(id),vec(cc[:,1]),nm,nm)

		id = ipts -mm*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,2]),nm,nm)

		id = ipts -(mm-1)*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,3]),nm,nm)

		id = ipts -1*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,4]),nm,nm)

		id = ipts + 1*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id), vec(cc[:,6]),nm,nm)

		id = ipts + (mm - 1)*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,7]),nm,nm)

		id = ipts + mm*ones(mm-2)
		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,8]),nm,nm)

		id = ipts + (mm + 1)*ones(mm-2)
    		AA = AA + sparse(vec(ipts),vec(id),vec(cc[:,9]),nm,nm)
    
		cc = nothing
	end
 
 

	B = 0*ones(nm)

	inn = Int(trunc(nn/2))
	iforce = mm*inn:mm:(mm*inn + 4*mm)
	B[iforce] = ones(size(iforce))

	bb = AA\B

	area = dr*(nn-3)*maximum(h)
	rrr = jLWavesInt(bb)
	BB = bb*sqrt(rrr*area)
	bref = BB[2*mm-1]
	fact = conj(bref)/abs(bref)
	BB = BB*fact

	@printf("wr, rrr = %.4e %.4eim  %.4e \n", real(w), imag(w), rrr)
	

	return rrr
end							# end of function







"""
		jLWavesCc(w)

	Compute coefficients for the coastal boundary for an island

K. Brink	4/1/2026
"""


function  jLWavesCc(w)

#   Compute coefficients for coastal boundary condition (avoiding very top and bottom grid points)
#	   Island version  
#
#   K. Brink 4/1/2026 based on Matlab of 9/2/2003
	
	global rl, m2, n2, f, vm, h, R, nn, mm, vmr, vmrr, icbc, dr, dt, hr, m2r, m2z, rrzero


	hh = h[1]

	th = -1 .+ dt*(0:(mm-3))

	wp = w .+ rl*vm[1,:]/rrzero
	wps = wp.*wp
	m2n2 = m2[1,:]./n2[1,:]
	f1 = f .+ 2*vm[1,:]/rrzero
	f2 = f .+ vmr[1,:] +vm[1,:]/rrzero
	f3 = f .+ 2*vm[1,:]/rrzero
	fs = f2 - m2n2.*m2[1,:]./f3
	G = f1.*fs - wp.*wp


	thx = -(hr[1]/hh)*th
	thz = 1/hh

	m2n2b = m2n2[1]
	f1b = f1[1]
	fsb = fs[1]
	wpb = wp[1]
	wpsb = wps[1]
	Gb = G[1]
	bfac = R[1]/Gb

	cc = 0*ones(mm-2,14)*(1 + im)

	drdt = dr/dt
	dtdr2 = 2*dt*dr
	dt2 = 2*dt
	drs = dr*dr



#       Closed boundary at r = r_zero

    	d1 = -im*Gb*hh*wp./G
    	d2 =   im*hh*Gb*wp.*m2n2./G
    	d3 = -im*rl*hh*Gb*f1./(rrzero*G)

    	d4 =  -bfac*(f1b*fsb + wpsb)
    	d5 = bfac*2*wpsb*m2n2b
    	d6 =  -bfac*(2*wpb*rl*f1b)/rrzero
     

    	C1 = d1
    	C2 = (d1.*thx + d2*thz)
    	C3 = d3
    
    	D4 = d4
    	D5 = d4*thx[1] + d5*thz
    	D6 = d6
    
    
    	cc[:,2] = -im*C1
    	cc[:,4] = -im*C2*drdt
    	cc[:,5] = im*C3*2*dr
    	cc[:,6] = im*C2*drdt
    	cc[:,8] = im*C1
    
    	cc[:,10] .= -im*D4
    	cc[:,11] .= -im*D5*drdt
    	cc[:,12] .= im*2*dr*D6
    	cc[:,13] .= im*D5*drdt
    	cc[:,14] .= im*D4
    
    
    
    #   Now scale it for consistency

     	cc  = im*cc*hh/(drdt*drdt)

	return cc
 end							# end of function





"""
		jLWavesOff(w)

	Compute coefficients for the offshore boundary condition
		for islands

K. Brink 4/1/2026
"""



function  jLWavesOff(w)

#   Compute the coefficients for the outer boundary condition
#   		for islands
#
#	K. Brink, 4/1/2026 based on Matlab of  1/1/2005

	global rl, m2, n2, f, vm, h, hr, nn, mm, vmr, vmrr, m2r, m2z, dr, dt, iobc, R, rrzero 

	hh = h[nn-2]

	rad = rrzero + dr*(nn-3)
	wp = w .+ rl*vm[nn-2,:]/rad
	wps = wp.*wp
	m2n2 = m2[nn-2,:]./n2[nn-2,:]
	m2n2r = m2r[nn-2,:]./n2[nn-2,:] - m2z[nn-2,:].*m2n2./n2[nn-2,:]
	f1 = f .+ 2*vm[nn-2,:]/rad
	f2 = (f1 .+ vmr[nn-2,:]) +vm[nn-2,:]/rad
	f3 = f .+ 2*vm[nn-2,:]/rad
	fs = f2 .- m2[nn-2,:].*m2n2./f3

	f1r = 2*(vmr[nn-2,:] - vm[nn-2,:]/rad)/rad
	f2r = vmrr[nn-2,:] + vmr[nn-2,:]/rad - vm[nn-2,:]/(rad^2)
	f3r = f1r
	fsr = f2r - 2*m2[nn-2,:].*m2r[nn-2,:]./(n2[nn-2,:].*f3)
	fsr = fsr + f3r.*m2[nn-2,:].*m2n2./(f3.*f3)
	fsr = fsr + m2[nn-2,:].*m2n2.*m2z[nn-2,:]./(n2[nn-2,:].*f3)
	wpr = rl*(vmr[nn-2,:] - vm[nn-2,:]/rad)/rad
	G = f1.*fs - wps
	Gr = f1r.*fs + f1.*fsr - 2*wp.*wpr



	thz = 1/hh

	drdt2 = 4*dr*dt
	dr2 = 2*dr
	drs = dr*dr

	cc = 0*ones(mm-2,14)*(1 + im)



	if iobc > 0.5
#   d(ru)/dr = 0 (Open) outer boundary condition 

 		b1 = -wp
		b2 = thz*wp.*m2n2
		b3 = wp.*Gr./G  -wpr -rl*f1/rad -wp/rad
		b4 = thz*(-(Gr./G).*wp.*m2n2 + wp.*m2n2./rad)
		b4 = b4 + thz*(wpr.*m2n2  + wp.*m2n2r)
		b5 = rl*f1.*Gr./(rad*G) - rl*f1r/rad


    		cc[:,1] = b2/drdt2
    		cc[:,2] =  b1/drs - b3/dr2
    		cc[:,3] = -cc[:,1]
    		cc[:,4] = -b4*0.5/dt
    		cc[:,5] = b5 - 2*b1/drs
    		cc[:,6] = -cc[:,4]
    		cc[:,7] = -cc[:,1]
    		cc[:,8] = b1/drs + b3/dr2
    		cc[:,9] = cc[:,1]

    
    
    		cc = cc*hh*hh*dt*dt
    
    
	else
#   closed outer BC: requires shallow water when R != 0

    		th = -1 .+ dt*(0:mm-3)
    		thr = -(hr[nn-2]/hh)*th
    		drdt = dr/dt

    		Gb = G[1]
    		wpb = wp[1]
   	 	f1b = f1[1]
   	 	fsb = fs[1]
   	 	RR = R[nn-2]
    		m2n2b = m2n2[1]
    
    		d1 = -im*hh*Gb*wp./G
   	 	d2 = im*hh*Gb*wp.*m2n2./G
    		d3 = -im*rl*Gb*hh*f1./(rad*G)
    
   	 	d4 = -RR*(f1b*fsb + wpb*wpb)/Gb
    		d5 = RR*(2*wpb*wpb*m2n2b)/Gb
    		d6 = -RR*(2*f1b*rl*wpb)/(rad*Gb)
    
    		C1 = d1
    		C2 = (d1.*thr + d2*thz)
    		C3 = d3
    
    		D4 = d4
    		D5 = d4*thr[1] + d5*thz
    		D6 = d6
    
    
    		cc[:,2] = -im*C1
    		cc[:,4] = -im*C2*drdt
    		cc[:,5] = im*C3*2*dr
    		cc[:,6] = im*C2*drdt
    		cc[:,8] = im*C1
    
    		cc[:,10] .= -im*D4
    		cc[:,11] .= -im*D5*drdt
    		cc[:,12] .= im*2*dr*D6
    		cc[:,13] .= im*D5*drdt
    		cc[:,14] .= im*D4

   
    #   Now scale it for consistency

     		cc  = cc*hh/(drdt*drdt)
	end
	return cc

end							# end of function






"""
	jLWavesBot(w)

	Computer coefficients for the bottom boundary condition
		ofr an island

K. Brink  4/1/2026  
"""


function jLWavesBot(w)

#   Compute coefficients for bottom boundary condition
#		for an island
#
# K. Brink 4/1/2026 

	global nn,  dr, dt, R, Rr, m2, m2r, m2z, n2, n2z, f, rl, vm, vmr, vmrr, h, hr, hrr, rrzero

	rad = rrzero .+ dr*(0:(nn-3))

	m2n2 = m2[:,1]./n2[:,1]
	wp = w .+ rl*vm[:,1]./rad
	wps = wp.*wp
	f1 = f .+ 2*vm[:,1]./rad
	f2 = f .+ vmr[:,1] + vm[:,1]./rad
	f3 = f1
	fs = f2 - m2n2.*m2[:,1]./f3
	G = f1.*fs - wp.*wp
	RG = R./G;
	Gp = f1.*fs + wp.*wp

	wpr = rl*(vmr[:,1] - vm[:,1]./rad)./rad
	m2n2r = m2r[:,1]./n2[:,1] - m2n2.*m2z[:,1]./n2[:,1]
	f1r = 2*(vmr[:,1] - vm[:,1]./rad)./rad
	f3r = f1r
	f2r = vmrr[:,1] + vmr[:,1]./rad - vm[:,1]./(rad.*rad)
	fsr = f2r - m2n2r.*m2[:,1]./f3 - m2n2.*m2r[:,1]./f3 + m2n2.*m2[:,1].*f3r./(f3.*f3)
	Gr = f1r.*fs + f1.*fsr - 2*wp.*wpr

	wpz = rl*m2[:,1]./(rad.*f3)
	m2n2z = m2z[:,1]./n2[:,1] - m2n2.*n2z[:,1]./n2[:,1]
	f1z = 2*m2[:,1]./(rad.*f3)
	f2z = (m2r[:,1] - m2[:,1].*f3r./f3 + m2[:,1]./rad)./f3
	f3z = f1z
	fsz = f2z -m2n2z.*m2[:,1]./f3 - m2n2.*m2z[:,1]./f3 + m2n2.*m2[:,1].*f3z./(f3.*f3)
	Gz = f1.*fsz + f1z.*fs - 2*wp.*wpz


	Grb = Gr - hr.*Gz
	m2n2rb = m2n2r - hr.*m2n2z
	f1rb = f1r - hr.*f1z
	f3rb = f1rb
	wprb = wpr - hr.*wpz
	fsrb = fsr - hr.*fsz
	GrG = Grb./G

	th = -1
	thz = ones(nn-2,1)./h
	thxz = -(hr./(h.*h))
	thx = (hr./h)
	thxbx = (hrr./h) - thx.*thx


	b1 = RG.*Gp
	b2 = -2*RG.*wps.*m2n2
	b3 = -im*wp.*m2n2 + im*wp.*hr +RG.*(-2*GrG.*Gp + wp.*wprb + f1.*fsrb + (Gp -2*rl*wp.*fs)./rad)
	b3 = b3 -(wp./G).*(wprb.*R + wp.*Rr) +(fs./G).*(f1rb.*R + f1.*Rr)
	b4 = 2*rl*f1.*wp.*RG./rad
	b5 = im*(wp./n2[:,1]).*G  + im*(wp.*m2n2).*(-hr + m2n2)
	b5 = b5 - (wp.*m2n2./(G.*f3)).*(f1.*(wprb.*R + wp.*Rr) +wp.*(f1rb.*R + f1.*Rr)) 
	b5 = b5 + RG.*(4*GrG.*wps.*m2n2 -wp.*(wprb.*m2n2 + wp.*m2n2rb))
	b5 = b5 + RG.*(-f1.*(2*wp.*wprb.*m2n2./f3 + wps.*m2n2rb./f3 - wps.*m2n2.*f3rb./(f3.*f3)))
	b5 = b5 + RG.*(m2n2.*(-2*wps + rl*wp.*fs + rl.*wp.*wps./f3)./rad)
	b6 = -(im*rl*f1./rad).*(m2n2 - hr) + (rl./(rad.*G)).*(f1.*(wprb.*R + Rr.*wp) + wp.*(f1rb.*R + Rr.*f1))
	b6 = b6 + RG.*(-4*rl*GrG.*f1.*wp +rl*wp.*f1rb + rl*f1.*wprb  -rl*rl*Gp./rad)./rad


	d1 = b1
	d2 = b1.*thx + b2.*thz
	d3 = b1.*thxbx + b3.*thx + b5.*thz + b2.*thxz
	d4 = b3 + b4
	d5 = b6

#   Scale for consistency

	d1 = dt*d1./(thz)
	d2 = dt*d2./(thz)
	d3 = dt*d3./(thz)
	d4 = dt*d4./(thz)
	d5 = dt*d5./(thz)

	drdt4 = 4*dr*dt
	drs = dr*dr
	dr2 = 2*dr
	dt2 = 2*dt

	cc= 0*ones(nn-2,9)*(1 + im)
	cc[:,1] = d2/drdt4
	cc[:,2] = d1/drs - d4/dr2
	cc[:,3] = -d2/drdt4
	cc[:,4] = -d3/dt2
	cc[:,5] = -2*d1/drs + d5
	cc[:,6] = d3/dt2
	cc[:,7] = -d2/drdt4
	cc[:,8] = d1/drs + d4/dr2
	cc[:,9] = d2/drdt4

	cc = im*cc

	return cc

	end								# end of function







"""
		jLWavesSbc(w)

	Compute coefficients for the surface boundary condition
		for an island

K.Brink  4/1/2026 
"""


function  jLWavesSbc(w)


#   Compute coefficients for use in the surface boundary condition
#		for an island
#   
#   
#	K. Brink 3/20/2026 

	global del, nn, mm, f, rl, vm, m2, n2, vmr, h, dr, dt, rrzero

	g = 980.

	cc = 0*ones(nn-2,5)*(1 + im)

	rad = rrzero .+ dr*(0:nn-3)

	wp = w .+ rl*vm[:,mm-2]./rad
	m2n2 = m2[:,mm-2]./n2[:,mm-2]
	f1 = f .+ 2*vm[:,mm-2]./rad
	f2 = f .+ vmr[:,mm-2] + vm[:,mm-2]./rad
	f3 = f .+ 2*vm[:,mm-2]./rad
	fs = f2 -m2[:,mm-2].*m2n2./f3
	G = f1.*fs -wp.*wp


	thz = ones(nn-2,1)./h
	thx = 0

	s1 = ((-wp./n2[:,mm-2]).*(G + m2[:,mm-2].*m2n2))./thz
	s2 = (wp.*m2n2)./thz
	s3 = (rl*f1.*m2n2./rad - (wp*del).*G/g)./thz

	c1 = s1.*thz + s2*thx

	cc[:,1]= -s2
	cc[:,2]= -c1*dr/dt
	cc[:,3]= 2*dr*s3
	cc[:,4] = c1*dr/dt
	cc[:,5] = s2

	cc = cc*dt/dr

	return cc

end						# end of function







"""
		jLWavesMatco(w,n)

	Compute coefficients for interior points for an island

K. Brink 4/1/26 
"""



function jLWavesMatco(w,n)

#  Get coefficients for the big array interior points
#		for an island
#  This works for a particular radial grid point n  
#
#  K. Brink 4/1/2026 based on Matlab of 5/18/2018

	global  mm, dt, dr, f, rl, vm, vmr, vmrr, m2, n2, n2z, m2z, m2r, h, hr, hrr, rrzero


	n = n-1
	rad =  rrzero + dr*(n - 1)

	wp = w .+ rl*vm[n,:]/rad

	wps = wp.*wp
	m2n2 = m2[n,:]./n2[n,:]
	m2n2r = m2r[n,:]./n2[n,:] - m2n2.*m2z[n,:]./n2[n,:]
	m2n2z = m2z[n,:]./n2[n,:] - m2n2.*n2z[n,:]./n2[n,:]
	f1 = f .+ 2*vm[n,:]/rad
	f2 = f .+ vmr[n,:] + vm[n,:]/rad 
	f3 = f .+ 2*vm[n,:]/rad
	f1r = 2*(vmr[n,:] - vm[n,:]/rad)/rad
	f1z = 2*m2[n,:]./(rad*f3)
	f2r = vmrr[n,:] + vmr[n,:]/rad - vm[n,:]/(rad*rad)
	f3r = 2*vmr[n,:]/rad  -2*vm[n,:]/(rad^2)	
	f2z = (m2r[n,:] - m2[n,:].*f3r./f3 + m2[n,:]./rad)./f3
	f3z = 2*m2[n,:]./(rad*f3)
	wpz = rl.*m2[n,:]./(rad*f3)
	wpr = rl*(vmr[n,:] - vm[n,:]/rad)/rad
	fs = f2 - (m2n2.*m2[n,:])./f3
	fsr = f2r -m2n2r.*m2[n,:]./f3 - m2n2.*m2r[n,:]./f3 + m2n2.*m2[n,:].*f3r./(f3.*f3)
	fsz = f2z -m2n2z.*m2[n,:]./f3 - m2n2.*m2z[n,:]./f3 + m2n2.*m2[n,:].*f3z./(f3.*f3)
	G = f1.*fs - wps
	Gr = f1.*fsr + fs.*f1r -2*wp.*wpr
	Gz = f1z.*fs + fsz.*f1 - 2*wp.*wpz



	Gamma = -Gr./G + (1/rad)*ones(mm-2) -m2n2z + m2n2.*Gz./G


	a1 = wp
	a2 = -2*wp.*m2n2
	a3 = wp.*(m2n2.*m2n2 + G./n2[n,:])
	a4 = wp.*Gamma + wpr + rl*f1/rad - m2n2.*wpz - rl*fs/rad
	a5 = -wp.*m2n2.*Gamma - rl*f1.*m2n2./rad + rl*wp.*wp.*m2n2./(rad*f3)
	a5 = a5 - wpr.*m2n2  - wp.*m2n2r
	a5 = a5 + m2n2.*(wpz.*m2n2  + wp.*m2n2z)
	a5 = a5 + G.*(wpz./n2[n,:] - wp.*n2z[n,:]./(n2[n,:].*n2[n,:]))
	a6 = (rl*Gamma.*f1/rad) + rl*(f1r/rad - f1/(rad*rad)) - rl*m2n2.*f1z/rad -rl*rl*wp/(rad*rad)



	th = -1 .+ dt*(0:(mm-3))
	thz = 1/h[n]
	thx = -hr[n]*th*thz
	thxx = -(hrr[n]/h[n])*th - 2*(hr[n]/h[n])*thx
	thxz = -hr[n]*thz*thz

	b1 = a1
	b2 = 2*thx.*a1 + a2.*thz
	b3 = (thx.*thx).*a1 + a3*thz*thz
	b3 = b3 + thz*thx.*a2
	b4 = a1.*thxx + a2.*thxz + a4.*thx + a5.*thz
	b5 = a4
	b6 = a6

	dt2 = 2*dt
	drdts2 = 2*dr/(dt*dt)
	drh = 2/dr
	drdt = dr/dt

	cc = 0*ones(mm-2,9)*(1 + im)
	cc[:,1] = b2/dt2
	cc[:,2] = b1*drh - b5
	cc[:,3] = -b2/dt2
	cc[:,4] = b3*drdts2 - b4*drdt
	cc[:,5] = -4*b1/dr + 2*dr*b6 - 2*b3*drdts2
	cc[:,6] = b3.*drdts2 + b4*drdt
	cc[:,7] = -b2/dt2
	cc[:,8] = b1*drh + b5
	cc[:,9] = b2/dt2

#   Scale for consistency
	cc = cc*h[n]*h[n]*dt*dt/dr

	return cc

end								# end of function







"""
		jLWavesInt(bb)

	Compute an (r, z)integral of the pressure squared
		for an island

K. Brink  4/1/2026 
"""


function  jLWavesInt(bb)

#   Compute the inverse response integral for the given iteration
#		for an island
#   bb comes in the form of an (nn*mm,1) complex array      5/5/03
#
# K. Brink 4/1/2026

	global nn, mm, dt, dr, h, rrzero

	r = rrzero .+ dr*(0:(nn-3))
	dz = h[1]*dt
	nlow = mm + 2
	nhigh = 2*mm - 1
	rr = 0.5*dz*r[1]*(abs(bb[nlow])^2 + abs(bb[nhigh])^2)
	rr = rr+ dz*r[1]*( 0.5*sum(abs.(bb[(nlow + 1):(nhigh - 1)]).^2))

	dz = h[nn-2]*dt
	ntop = (nn-1)*mm - 1
	nbot = (nn-2)*mm + 2
	rt = 0.5*dz*r[nn-2]*(abs(bb[ntop])^2 + abs(bb[nbot])^2)
	rt = rt + dz*r[nn-2]*( 0.5*sum(abs.(bb[(nbot+1):(ntop-1)]).^2))

	rrr = 0

	for n = 3:(nn-2)
	    dz = h[n-1]*dt      
            ntop = n*mm - 1
	    nbot = (n-1)*mm + 2
	    rs = 0.5*dz*r[n-1]*(abs(bb[nbot])^2 + abs(bb[ntop])^2)
	    rs = rs + dz*r[n-1]*(sum(abs.(bb[(nbot+1):(ntop-1)]).^2))
	    rrr = rrr + rs
	end

	rrr = dr*(rrr + rt + rr)
	rrr = 1/rrr

	return rrr

	end						# end of function







"""
		jLWavesEngDiag(war,BB)

	Compute energy diagnostics for islands

K. Brink	4/1/2026 
"""

	


function jLWavesEngDiag(war,BB)

#   Compute energy diagnostics: most useful for unstable waves
#
#   	island version 
#
# K. Brink 4/1/2026 

	global nn, mm, dr, dt, f, rl, del, R, h, n2, m2, vm, vmr,rrzero 

	con1 = 0
	con2 = 0
	con3 = 0
	con4 = 0
	con5 = 0

	epes = 0*(1 + im)
	eke = 0*copy(epes)
	epe = 0*copy(epes)

	
	g = 980
	rhob = 1.03
	gr = g*rhob



	utup = jLWavesUvwr(war,BB)
	u = utup[1]
	v = utup[2]
	wvel = utup[3]
	rho = utup[4]

	uc = conj(u)
	vc = conj(v)
	wvelc = conj(wvel)
	rhoc = conj(rho)

	w = war
	m2n2 = m2./n2
	dr2 = 2*dr


	for n = 1:nn-2
    		rad = rrzero + dr*(n - 1)
    		dz = dt*h[n]
    		drl = dr
		if n == nn-2 || n == 1
    			drl = dr/2
		end		
    
    		itop = ((n+1)*mm) - 1
    		epes = epes + rad*drl*BB[itop]*conj(BB[itop])
    
    		xx = u[n,:].*uc[n,:] + v[n,:].*vc[n,:]
    		yy = (rho[n,:].*rhoc[n,:])./n2[n,:]
    		eke = eke + rad*drl*dz*(sum(xx) - (xx[1] +xx[mm-2])/2)
    		epe = epe + rad*drl*dz*(sum(yy) - (yy[1] +yy[mm-2])/2)
    
    		s1 = (vc[n,:].*u[n,:] + uc[n,:].*v[n,:]).*vmr[n,:]
    		s2 = (vc[n,:].*wvel[n,:] + (wvelc[n,:].*v[n,:])).*m2[n,:]/f
    		s3 = wvelc[n,:].*rho[n,:] + wvel[n,:].*rhoc[n,:]
    		s5 = (u[n,:].*rhoc[n,:] + uc[n,:].*rho[n,:]).*m2n2[n,:]
    
    		con1 = con1 + rad*drl*dz*(sum(s1) - (s1[1] +s1[mm-2])/2)
    		con2 = con2 + rad*drl*dz*(sum(s2) - (s2[1] +s2[mm-2])/2)
    		con3 = con3 + rad*drl*dz*(sum(s3) - (s3[1] +s3[mm-2])/2)
   		con5 = con5 + rad*drl*dz*(sum(s5) - (s5[1] +s5[mm-2])/2)
    
   		 wp = w .+ rl*vm[n,1]/rad
   		 wps = wp*wp
    		f1 = f + 2*vm[n,1]/rad
    		f2 = f + vmr[n,1] + vm[n,1]/rad
    		f3 = f + 2*vm[n,1]/rad
    		fs = f2 - m2[n,1].*m2n2[n,]./f3
    		G = f1.*fs - wps
    
    		rG = R[n]/G
    
     		imain = n*mm+2
     		pr = (BB[imain+mm] - BB[imain-mm])/dr2
   
     		pbr = pr
     		pbrc = conj(pbr)
    
    		 Ue = rG*(-im*wp*u[n,1] - f1*v[n,1])
    		 Uec = conj(Ue)
     		Ve = rG*(fs*u[n,1] - im*wp*v[n,1])
     		con4 = con4 + rad*drl*(pbr*Uec + pbrc*Ue)
    
 	end

	eke = real(eke*rhob/2)
	epes = real(epes*del/(gr*2))
	epe = real(0.5*epe*g*g/rhob)

	con1 = real(-con1/2)
	con2 = real(-con2/2)
	con3 = real(-con3*g/2)
	con4 = real(con4/2)
	con5 = real(con5*g/2)

	ss = sum(con1 + con2 + con5 + con4)
	if ss == 0
	    con3 = 0
	end

#   Plot results



	fig10 = Figure()
	

	#colsize!(fig6.layout, 1, Aspect(1, 1.5))
				
	ax1 = Axis(fig10[1,1])	


#				 Define arrows
	arx = [-3, 0, -3]
	ary = [-3, 0, 3]
	arvx = [-2, 0, 2]
	arvy = [-3, 0, -3]

	x = [35, 75, 75, 35, 35]
	y = [60, 60, 90, 90, 60]

	lines!(ax1,x,y, color = :black)
	limits!(ax1,0, 100, 0, 100)

	y = y .- 50
	lines!(ax1,x,y, color = :black)
	text!(ax1, 50, 80, text = "EKE", color = :black, fontsize = 20)
	stt = string(@sprintf("%.3e",eke))
	text!(45,70, text = stt, color = :black,fontsize = 20)
	text!(50,30,text = "EPE", color = :black, fontsize = 20)
	eee = epe + epes
	stt = string(@sprintf("%.3e",eee ))
	text!(45,19, text = stt, color = :black,fontsize = 20)

	x = [15, 35]
	y = [85, 85]
	lines!(ax1,x,y, color = :black)
	if con1 >= 0
 		lines!(ax1,arx .+ 35,ary .+ 85, color = :black)
	else
		lines!(ax1, -arx .+ 15, ary .+ 85, color = :black)
	end
	text!(ax1, 1,90,text = "Barotropic MKE to EKE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con1))
	text!(ax1,19,80,text = stt, color = :black, fontsize = 15)

	y = y .- 20
	lines!(ax1,x,y, color = :black)
	if con2 >= 0
 		lines!(ax1,arx .+ 35,ary .+ 65, color = :black)
	else
		lines!(ax1, -arx .+ 15, ary .+ 65, color = :black)
	end
	text!(ax1, 5,70,text = "Shear MKE to EKE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con2))
	text!(ax1,17,60,text = stt, color = :black, fontsize = 15)

	y = [25, 25]
	lines!(ax1,x,y, color = :black)
	if con5 >= 0
 		lines!(ax1,arx .+ 35,ary .+ 25, color = :black)
	else
		lines!(ax1, -arx .+ 15, ary .+ 25, color = :black)
	end
	text!(ax1, 13, 30,text = "MPE to EPE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con5))
	text!(ax1,19,20,text = stt, color = :black, fontsize = 15)

	x = [55, 55]
	y = [40, 60]
	lines!(ax1,x,y, color = :black)
	if con3 >= 0
 		lines!(ax1,arvx .+ 55,arvy .+ 60, color = :black)
	else
		lines!(ax1, arvx .+ 55, -arvy .+ 40, color = :black)
	end
	text!(ax1, 34,50,text = "EPE to EKE", color = :black, fontsize = 18)
	stt = string(@sprintf("%.3e",con3))
	text!(ax1,58,50,text = stt, color = :black, fontsize = 15)

	x = [75, 95]
	y = [75, 75]
	lines!(ax1,x,y, color = :black)
	if con4 <= 0
 		lines!(ax1,arx .+ 95,ary .+ 75, color = :black)
	else
		lines!(ax1, -arx .+ 75, -ary .+ 75, color = :black)
	end
	text!(ax1, 80,80,text = "Dissipation", color = :black, fontsize = 20)
	stt = string(@sprintf("%.3e",con4))
	text!(ax1,78,70,text = stt, color = :black, fontsize = 15)
	#	This should be a negative number



	stt = string(@sprintf("%.3e %.3eim",real(war), imag(war)))
	ax1.title = string("Energy diagnostics for ω = ", stt, " 1/sec")
		
	#hidedecorations!(ax1, label = false, ticklabels = false, ticks = false )
	hidedecorations!(ax1)

	display(fig10)


end							# end of function







"""
		jLWavesUvwr(ww,BB)

	Compute u, v, w, rho once the wave mode is known
		for an island

K. Brink 4/1/2026 

"""



function jLWavesUvwr(ww,BB)

#   Compute u, v, w and rho once the pressure is known
#		for an island
#   ww is the frequency and BB is the pressure vector
#
# K. Brink 4/1/2026 


	global nn, mm, dr, dt, rl, f, vm, vmr, n2, m2, h, hr, rrzero

	w = ww

	u = 0*ones(nn-2,mm-2)*(1 + im)
	v = copy(u)
	wvel = copy(u)
	rho = copy(u)
	pp = copy(u)
	ptop = 0*ones(nn-2)*(1 + im)
	pbot = copy(ptop)

# Max ratio of real/imag or imag/real to allow plotting of the smaller 
	ratnplot = 1000

	g = 980
	rhob = 1.03
	dt2 = 2*dt
	dr2 = 2*dr
	th = -1 .+ dt*(0:(mm-3))


	for n = 1:(nn-2)
	    rad = rrzero + dr*(n - 1)
	    irow = (n*mm + 2):((n+1)*mm -1)
	    wp = w .+ rl*vm[n,:]/rad
	    wps = wp.*wp
	    m2n2 = m2[n,:]./n2[n,:]
	    f1 = f .+ 2*vm[n,:]/rad
	    f2 = f .+ vmr[n,:] + vm[n,:]/rad
	    f3 = f .+ 2*vm[n,:]/rad
    
	    fs = f2 - m2[n,:].*m2n2./f3
	    G = f1.*fs - wp.*wp
	    thx = -(hr[n]/h[n])*th
    
	    ptop[n] = BB[(n+1)*mm-1]
	    pbot[n] =  BB[(n*mm +2)]
	    p = BB[irow]
	    pt = (BB[irow .+ 1] - BB[irow .- 1])/dt2
	    prp = (BB[irow .+ mm] - BB[irow .- mm])/dr2
       
	    pr = prp + thx.*pt
	    pz = pt/h[n]
    

	    u[n,:] = (-im*wp.*pr - (im*rl*f1./rad).*p + im*wp.*m2n2.*pz)./G
	    v[n,:] = (fs.*pr + (rl*wp/rad).*p - (wps.*m2n2./f3).*pz)./G
   
	    wvel[n,:] = (((-im*wp./n2[n,:])).*pz) - m2n2.*u[n,:]
	    rho[n,:] = -pz/g

	    pp[n,:] = p
    
	end

	u = u/rhob
	v = v/rhob
	wvel = wvel/rhob

	utup = (u, v, wvel, rho)

#  Plot results

	rgr = ones(nn-2,mm-2)
	zgr = copy(rgr)
    
 	for n = 1:nn-2       
        	rtemp =  (rrzero .+ dr*(n - 1))/1e5 
        	rgr[n,:]= rtemp*ones(1,mm-2)
        	zgr[n,:] = -h[n]*(1 .- dt*(0:mm-3))/100
 	end



	hh = -h/100
	rpl = (rrzero .+ dr*(0:nn-3))/1e5
	ax = [0 maximum(rpl) minimum(hh) 0]
	z = th*maximum(h)/100
	xf = [0 rpl' maximum(rpl) 0]
	zf = [minimum(hh) hh' minimum(hh) minimum(hh)]
	rdep = vec(rgr[:,1])
	rmax = maximum(rdep)
	rmin = minimum(rdep)
	hmax = maximum(h)


	fig5 = Figure()
	

	colsize!(fig5.layout, 1, Aspect(1, 1.5))
				
	ax1 = Axis(fig5[1,1], title = "uplot")				# top panel
	ax3 = Axis(fig5[2,1], title = "wplot")				# lower
	ax2 = Axis(fig5[1,2], title = "vplot")
	ax4 = Axis(fig5[2,2], title = "rplot")


	rmm = maximum(maximum(abs.(real(pp))))
	amm = maximum(maximum(abs.(imag(pp))))
	rata = amm/rmm
	pmult = 1.
	if rata < (1/ratnplot)
		pmult = NaN
	end


#		plot u
	rmm = maximum(maximum(abs.(real(u))))
	amm = maximum(maximum(abs.(imag(u))))
	rata = rmm/amm
	umult = 1.
	if rata < (1/ratnplot)
		umult = NaN
	end

	civ = jHWavesCCon(u,10)
	
	zmax = maximum(-hh)
	yupper = hh
	nzz = length(yupper)
	ylower = 0*yupper -zmax*ones(nzz)
	band!(ax1,rdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax1,rdep,hh, color = (:blue,0.7))
	xlims!(ax1,rmin,rmax)
	ylims!(ax1,-hmax/100,0)
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)
 	ax1.ylabel = "z [m]"
	ax1.xlabel = "r [km]"
	ax1.title = "u"

	contour!(ax1,rgr,zgr,umult*real(u), levels = civ, color = :black)
	contour!(ax1,rgr,zgr,umult*real(u),levels = vec([0 0]), linewidth = 2, color = :black)
	contour!(ax1,rgr,zgr,imag(u), levels = civ, color = :red)
	contour!(ax1,rgr,zgr,imag(u),levels = vec([0 0]), linewidth = 2, color = :red)

	



#		plot v
	civ = jHWavesCCon(v,10)

	
	rmm = maximum(maximum(abs.(real(v))))
	amm = maximum(maximum(abs.(imag(v))))
	rata = amm/rmm
	vmult = 1.
	if rata < (1/ratnplot)
		vmult = NaN
	end


	band!(ax2,rdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax2,rdep,hh, color = (:blue,0.7))
	xlims!(ax2,rmin,rmax)
	ylims!(ax2,-hmax/100,0)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)
 	ax2.ylabel = "z [m]"
	ax2.xlabel = "r [km]"
	ax2.title = "v"

	contour!(ax2,rgr,zgr,real(v), levels = civ, color = :black)
	contour!(ax2,rgr,zgr,real(v),levels = vec([0 0]), linewidth = 2, color = :black)
	contour!(ax2,rgr,zgr,vmult*imag(v), levels = civ, color = :red)
	contour!(ax2,rgr,zgr,vmult*imag(v),levels = vec([0 0]), linewidth = 2, color = :red)


#		plot w

	civ = jHWavesCCon(wvel,10)
	

	rmm = maximum(maximum(abs.(real(w))))
	amm = maximum(maximum(abs.(imag(w))))
	rata = rmm/amm
	wmult = 1.
	if rata < (1/ratnplot)
		wmult = NaN
	end

	band!(ax3,rdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax3,rdep,hh, color = (:blue,0.7))
	xlims!(ax3,rmin,rmax)
	ylims!(ax3,-hmax/100,0)
	hidedecorations!(ax3, grid = true, ticklabels = false, label = false)
 	ax3.ylabel = "z [m]"
	ax3.xlabel = "r [km]"
	ax3.title = "w"

	contour!(ax3,rgr,zgr,imag(wvel), levels = civ, color = :red)
	contour!(ax3,rgr,zgr,imag(wvel),levels = vec([0 0]), linewidth = 2, color = :red)
	contour!(ax3,rgr,zgr,wmult*real(wvel), levels = civ, color = :black)
	contour!(ax3,rgr,zgr,wmult*real(wvel),levels = vec([0 0]), linewidth = 2, color = :black)



#		plot rho
	civ = jHWavesCCon(rho,10)

	rmm = maximum(maximum(abs.(real(rho))))
	amm = maximum(maximum(abs.(imag(rho))))
	rata = amm/rmm
	rmult = 1.
	if rata < (1/ratnplot)
		rmult = NaN
	end


	band!(ax4,rdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax4,rdep,hh, color = (:blue,0.7))
	xlims!(ax4,rmin,rmax)
	ylims!(ax4,-hmax/100,0)
	hidedecorations!(ax4, grid = true, ticklabels = false, label = false)
 	ax4.ylabel = "z [m]"
	ax4.xlabel = "r [km]"
	ax4.title = "rho"

	contour!(ax4,rgr,zgr,real(rho), levels = civ, color = :black)
	contour!(ax4,rgr,zgr,real(rho),levels = vec([0 0]), linewidth = 2, color = :black)
	contour!(ax4,rgr,zgr,rmult*imag(rho), levels = civ, color = :red)
	contour!(ax4,rgr,zgr,rmult*imag(rho),levels = vec([0 0]), linewidth = 2, color = :red)



	display(GLMakie.Screen(),fig5)


#       Plot surface values

#		Scale so that abs(v(0,0) = 1
	ssff = maximum(abs.(v[:,end]))
	ratp = 1/ssff

	fig6 = Figure()
	

	colsize!(fig6.layout, 1, Aspect(1, 1.5))
				
	ax1 = Axis(fig6[1,1], title = "pplot")				# top panel
	ax3 = Axis(fig6[3,1], title = "uplot")				# lower
	ax2 = Axis(fig6[2,1], title = "vplot")
	ax4 = Axis(fig6[4,1], title = "hplot")




#		plot surface p
	lines!(ax1,vec(rdep),real(vec(ptop*ratp)), color = :black)
	lines!(ax1,vec(rdep),pmult*imag(vec(ptop*ratp)), color = :red)
	lines!(ax1,vec(rdep),vec(0*rdep),color = :black)
	xlims!(ax1,rmin,rmax)
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)
 	ax1.ylabel = "p [dyne/cm^2]"
	ax1.xlabel = "r [km]"
	ax1.title = "Pressure (all amplitudes arbitrary, but consistent)"


#		plot v
	lines!(ax2,vec(rdep),real(vec(v[:,end]*ratp)), color = :black)
	lines!(ax2,vec(rdep),vmult*imag(vec(v[:,end]*ratp)), color = :red)
	lines!(ax2,vec(rdep),vec(0*rdep), color = :black)
	xlims!(ax2,rmin,rmax)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)
 	ax2.ylabel = "v [cm/s]"
	ax2.xlabel = "r [km]"
	ax2.title = "Surface v"


#		plot u
	lines!(ax3,rdep,umult*real(u[:,end]*ratp), color = :black)
	lines!(ax3,rdep,imag(u[:,end]*ratp), color = :red)
	lines!(ax3,rdep,0*rdep,color = :black)	
	xlims!(ax3,rmin,rmax)
	hidedecorations!(ax3, grid = true, ticklabels = false, label = false)
 	ax3.ylabel = "u [cm/s]"
	ax3.xlabel = "r [km]"
	ax3.title = "Surface u"

#		plot depth
	
	band!(ax4,rdep,ylower,yupper,color = (:blue,0.7))
	lines!(ax4,rdep,hh, color = (:blue,0.7))
	xlims!(ax4,rmin,rmax)
	ylims!(ax4,-hmax/100,0)
	
	tt = maximum(maximum(abs.(vm)))
	if tt > 0.1
		contour!(ax4,rgr,zgr,vm, levels = 8, color = :black)
	end

	hidedecorations!(ax4, grid = true, ticklabels = false, label = false)
 	ax4.ylabel = "z [m]"
	ax4.xlabel = "r [km]"
	ax4.title = "Depth profile and mean velocity"
	
	display(GLMakie.Screen(),fig6)

	return utup




end						# end of function








"""
	ci = jHWavesCCon(u,nc)

	Define contour intervals, given the range
	Done so that the real and imaginary parts have the same interval

K. Brink
"""

function jHWavesCCon(u, nc)
	


"""	
	jHWavesCCon(u, nc)

	function to make consistent contour intervals for a complex array u
		for use with jHWavesC and other codes with complex frequencies
	nc = # of contours
	Returns an array of contour intervals
K. Brink 3/15/2026
"""

	urma = maximum(maximum(real(u)))
	urmi = minimum(minimum(real(u)))
	uima = maximum(maximum(imag(u)))
	uimi = minimum(minimum(imag(u)))
	umi = uimi
	if urmi < uimi
		umi = urmi
	end
	uma = uima
	if urma > uima
		uma = urma
	end
	if umi > 0
		umi = 0
	end	
	if uma < 0
		uma = 0
	end
	ci = (uma-umi)/nc

	cc = umi:ci:uma
	(zz,ind) = findmin(abs.(cc))
	cont = cc .- zz
		
	return cont
end








"""
	jLWavesNorm(pin)

	Normalize an island mode

K. Brink 4/1/2026 
"""

	function jLWavesNorm(pin)

#	Normalize an island wave
#		assumes that the far boundary is open
#
#  K. Brink  4/1/2026 based on Matlab of 6/6/2018


	global nn, mm, dr, dt, h, hr, BB, rrzero



	dz = dt*h[1]
	pc = pin[1,:]
	con1 = sum(abs2.(pc)) - 0.5*(abs2(pc[1]) + abs2(pc[mm-2]))
	con1 = con1*dz

	pb = hr.*abs2.(pin[:,1])
	con2 = sum(pb) - 0.5*(pb[1] + pb[end])
	con2 = con2*dr

	po = pin[end,:]
	con3 = sum(abs2.(po)) - 0.5*(abs2(po[1]) + abs2(po[mm-2]))
	con3 = con3*dt*h[end]
	if iobc > 0.5
		con3 = 0
	end

	con = sqrt(abs(con1 + con2 - con3))

	p =  pin/con
	BB = BB/con

	return p

	end							# end of function




"""
	jLWavesSetup()

	Create the input tuple for use with jLWavesM (islands)

K.H. Brink 4/1/2026
"""


function jLWavesSetup()

#   Create an input tuple for "jLWavesM" (islands), given inputs
#   Call as
#       arr = jLWavesSetup()
#   In response to queries, you can enter either numbers or array names
#
#   K.H. Brink  4/1/2026

	println(' ')
	println("      Island-Trapped waves with complex frequency")
	println(' ')
	println("This function will ask you a sequence of questions ")
	println( "     that are used to build the input tuple.     ")
	println(' ')
    

	nn = inputi("How many total gridpoints do you want in the radial direction? (nn) ")
	if nn <5
    		println(' ')
    		println("Error! Must have nn > 4")
		return NaN
	end

	mm = inputi("How many total gridpoints do you want in the vertical? (mm) ")
	if mm < 5
    		println(' ')
    		println("Error! Must have mm >4")
		return NaN
	end
	arr = (nn, mm)

	wg = inputf("Enter first guess at real and complex parts of frequency (2-elemnent array) (wg) (1/sec) ")
		if length(wg) != 2
			println("Error: frequency should be an array of length 2!")
			return NaN
		end
	del = inputi("Enter 0 for a rigid lid, 1 for a free surface (del) ")
	icbc = NaN
	iobc = inputi("Enter 0 for a closed r =rmax boundary, 1 for an open (iobc) ")
	ilw = NaN
	f = inputf("Enter the Coriolis parameter (f) (rad/sec) ")
	rmax = inputf("Enter the maximum radius (rmax) (km) ")
	eps = inputf("Enter the nominal fractional accuracy for the frequency solution (eps) ")
	npts = inputi("Enter the number of frequencies to be computed (npts) ")

	
	
	arr = (arr..., wg, del, icbc, iobc, ilw, f, rmax, eps, npts)					# 11 items



	nwz = inputi("Enter the first azimuthal wavenumber to use (dnw) (integer) ")
	if npts < 1.1
    		dnw = nwz
	else
    		dnw = inputi("Enter the wavenumber increment to use after nwz (dnw) (integer) ")
	end
	arr = (arr..., nwz,dnw)										# 13 items

#   Read in depth
	println(' ')
	println("		When entering topography, the first radius falls at the coast")
	println(' ')

	ndep = inputi("How many radius, depth pairs will you provide (ndep >=1) ")
	
	rdep = inputf("Array of radial distances for depth values (rdep in km) (dimension ndep) ")
	if ndep != length(rdep)
		println(' ')
		println("Error! array length must equal ndep!")
		return NaN
	end
	depr = inputf("Array of depths corresponding to rdep (depr in m) ")
	if ndep != length(depr)
		println(' ')
		println("Error! array length must equal ndep!")
		return NaN
	end
	if rdep[1] < 0.000001
		println(' ')
		println("Error! The island radius must be > 0!")
		return NaN
	end
	if rdep[end] < rmax
		ndep = ndep + 1
		xdep = [rdep rmax]
		depr = [depr depr[end]]
	end

	if ndep ==1
		ndep = 2
		rdep = [rdep rmax]
		depr = [depr depr]
	end
	iqq = findall(depr -> depr < 0, depr)
	if length(iqq) > 0.5
		println(' ')
		println("Error! All depths must be > 0!")
		return NaN
	end
	arr = (arr..., ndep, rdep, depr)							# 16 items
		


#       Read in bottom friction
	nR = inputi("Number of radius, bottom friction pairs to read (nR) ")
	if nR == 0
		zc = 0
		nR = 2
		rR = [0 rmax]
		Rread = [0. 0.]
	else
    		rR = inputf("Radial distances for bottom friction values (rR in km) ")
    		Rread = inputf("Array of bottom friction values corresponding to rR (rr in cm/sec) ")
	        if length(rR) != nR
        		println(' ')
        		println("Error! Array size must match nR!")
        		return NaN
    		end
    		if length(Rread) != nR
        		println(' ')
        		println("Error! Array size must match nR!")
        		return NaN
	    	end
		if rR[1] > rdep[1]
			nR = nR + 1
			rR = [rdep[1] rR]
			Rread = [Rread[1] Rread]
		end
		if rR[end] < rmax
			nR = nR + 1
			rR = [rR rmax]
			Rread = [Rread Rread[end]]
		end
								# 19 items
	end
    	arr = (arr..., nR, rR, Rread)	
  

#       Read in base-state stratification
	zr = inputf("Depth increment for input of Nsquared values? (zr in m) ")
	alph = inputf("Exponential tail length for Nsquared extrapolation (alph in km) ")
	nsqr = inputf("Nsquared values starting at the surface (nsqr in rad^2/sec^2) (nnsq values) ")
	nnsq = length(nsqr)

	arr = (arr..., nnsq, zr, alph, nsqr)							# 23 items
	
#       Read in mean flow

	vzero = inputf("Input peak value of mean azimuthal flow (vzero: cm/sec) ")


	if abs(vzero) > 0.001
    		rzero = inputf(" Input radial distance to peak mean flow (km) ")
    		zzero = inputf(" Input depth of peak mean flow (m) ")
    		zscaled = inputf(" Downward exponential scale of mean flow? (m) ")
    		zscaleup = inputf(" Upward exponential scale of mean flow? (m) ")
    		xscaleoff = inputf(" Offshore exponential scale of mean flow? (km) ")
    		xscaleon = inputf(" Onshore exponential scale of mean flow? (km) ")
    		kk = inputi(" Enter 1 for undisturbed Nsquared offshore, 0 for onshore ")
		nv = 1
	else
		rzero = 20.
		zzero = 0.
		zscaled = 10000.
		zscaleup = 10000.
		xscaleoff = 1000.
		xscaleon = 1000.
		kk = 1
		nv = 0
	end

    	arr = (arr..., vzero, rzero, zzero, zscaled, zscaleup, xscaleoff, xscaleon, kk)		# 31 items

#	Miscellaneous
	ipause = inputi("Enter 0 to skip pauses to see graphics, 1 to see graphics during execution ")
	arr = (arr..., ipause)									# 32 items total



#	This completes the reading in of information
#		Now for some graphics

	
	fig = Figure()
	display(GLMakie.Screen(),fig)

	colsize!(fig.layout, 1, Aspect(1, 1.8))
				
	ax1 = Axis(fig[1,1], title = "rplot")				# top panel
	ax2 = Axis(fig[2,1], title = "vplot")				# lower

	lines!(ax1,vec(rR),vec(Rread))
	zc = maximum(Rread)
	xlims!(ax1,rdep[1],rmax)			
	if zc < 1e-8
		ylims!(ax1,0,0.1)

		xtex = rdep[1] + (rmax-rdep[1])/3
		text!(ax1,xtex,0.03, text = "R = 0")
	else
		ymaxx = 1.1*maximum(Rread)
		ylims!(ax1,0, ymaxx)
	end
	hidedecorations!(ax1, grid = true, ticklabels = false, label = false)
	ax1.ylabel = "R [cm/s]"
	ax1.title = "Friction Coefficient"

	zmax = maximum(depr)
	yupper = -depr
	nzz = length(yupper)
	ylower =  -zmax*ones(nzz)
	band!(ax2,vec(rdep),vec(ylower),vec(yupper),color = (:blue,0.7))
	lines!(ax2,vec(rdep),vec(-depr), color = (:blue,0.7))
	xlims!(ax2,rdep[1],rmax)
	ylims!(ax2,-zmax,0)
	hidedecorations!(ax2, grid = true, ticklabels = false, label = false)
 	ax2.ylabel = "z [m]"
	ax2.xlabel = "r [km]"
	vstr = @sprintf("Depth and Mean v (vzero = %.1f cm/sec)", vzero )
	ax2.title = vstr

	if nv == 0
		xtex = rdep[1] + 2*(rmax-rdep[1])/3		
		text!(ax2,xtex,-0.5*zmax, text = "Mean v = 0")
	else
		nnn = 41
	    	mmm = 31
	        r = 0:(rmax/(nnn-1)):rmax
    	        z = -zmax:(zmax/(mmm-1)):0
    		zzero = -zzero

    		xfactoff = xscaleoff^2
    		xfacton = xscaleon^2
    		zfactu = zscaleup^2
    		zfactd = zscaled^2
		
    
    		vtemp = NaN*ones(mmm,nnn)
    
    		ii = findall(r-> r < rzero,r)
    		iic = findall(r -> r >=  rzero,r)
    		jj = findall(z -> z < zzero,z)
    		jjc = findall(z -> z >= zzero,z)
 
    		if isempty(size(ii)) != 1
        	   for n = ii[1]:maximum(ii)
			
            		if isempty(size(jj)) != 1
				sjj = length(jj)
				zzz = zzero*ones(sjj)
               			 vtemp[jj,n] = (exp(-((r[n]-rzero).*(r[n]-rzero)/xfacton)))*(exp.(-(z[jj]-zzz).*(z[jj]-zzz)/zfactd))
            		end
           	 	if isempty(size(jjc)) != 1
				sjj = length(jjc)
				zzz = zzero*ones(sjj)
                		vtemp[jjc,n] = (exp(-((r[n]-rzero).*(r[n]-rzero)/xfacton)))*(exp.(-(z[jjc]-zzz).*(z[jjc]-zzz)/zfactu))
            		end
			
        	   end
    		end
    
    		if isempty(size(iic)) != 1
        	   for n = minimum(iic):nnn
            		if isempty(size(jj)) != 1
				sjj = length(jj)
				zzz = zzero*ones(sjj)
                 		vtemp[jj,n] = (exp(-((r[n]-rzero).*(r[n]-rzero)/xfactoff)))*(exp.(-(z[jj]-zzz).*(z[jj]-zzz)/zfactd))
            		end
            		if isempty(size(jjc)) != 1
				sjj = length(jjc)
				zzz = zzero*ones(sjj)
                 		vtemp[jjc,n] = (exp(-((r[n]-rzero).*(r[n]-rzero)/xfactoff)))*(exp.(-(z[jjc]-zzz).*(z[jjc]-zzz)/zfactu))
            		end
        	   end
   	 	end
		vtemp = vzero*vtemp
    
    
    		contour!(ax2,r, z, vtemp', levels = 5, color = :black, labels = false)
    		band!(ax2,vec(rdep),vec(ylower),vec(yupper),color = (:blue,0.7))
	end
	display(fig)

	vtemp = nothing
	return arr
end									# end of function



"""
	jLWavesFinch(arr)

	Make changes to input tuple for complex frequency
	Returns a revised tuple
	For island-trapped waves

K.H. Brink 4/2/2026
"""


function  jLWavesFinch(arr)

#   Function to modify an existing jLWaves (island) input tuple conveniently
#
#   K.H. Brink,  4/2/2026

	println(' ')
	println("     Function to create modified input tuple for an island")
	println(" ")

	xx = """First you need to select what you want to change.
                                                      
          Options are:                                
            g:  Grid size                     
            w:  Initial Frequency guess                                 
            b:  Boundary conditions
            f:  Coriolis parameter                     
            x:  Domain size                        
            d:  Dispersion curve definition
            e:  Nominal accuracy               
            h:  Depth profile                          
            R:  Bottom friction
            n:  Stratification
            v:  Mean flow field
            p:  Pauses during execution
                                                      
     Any arrays are row arrays, not column arrays.  

     Please type in the appropriate letter and hit return  """  
                                                      
 	println(xx)

	stinn = readline()
	stin = stinn[1]
	println(' ')   

	if stin == 'g'
    		println("Old [nn   mm] = ", arr[1], ",  ", arr[2])
    		iii = inputi(" Enter new [nn  mm] ")
		nn = iii[1]
		mm = iii[2]
    		narr = ( nn, mm, arr[3:end]...)

	elseif stin == 'w'
    		println("Old frequency guess (real and complex parts) (sec^-1) ",  arr[3])
    		www = inputf(" Enter new real and complex frequency parts (2-element array) (sec^-1) ")
		if length(www) != 2
			println("Error: frequency should be an array of length 2!")
			return
		end
		narr = (arr[1:2]..., www, arr[4:end]...)
    

	elseif stin == 'b'
    		println("Old boundary conditions")
    		if arr[4] == 0
        		println("   Rigid lid")
    		else
        		println("   Free surface")
    		end

    		if arr[6] == 0
        		println("   Closed boundary at r = rmax")
    		else
        		println("   Open boundary at r = rmax")
    		end
    		ar1 = inputi(" Enter 0 for rigid lid, 1 for free surface ")
   		ar2 = NaN
    		ar3 = inputi(" Enter 0 for closed offshore boundary, 1 for open ")
    		narr = (arr[1:3]..., ar1, ar2, ar3, arr[7:end]...)

		if ar3 > 0.5
			println(' ')
			println("Note: bottom must be flat at r = rmax")
		end

	elseif stin == 'f'
    		println("Old Coriolis parameter = ", arr[8])
    		www = inputf(" Enter new Coriolis parameter  (sec^-1) ")
		narr = (arr[1:7]..., www, arr[9:end]...)

	elseif stin == 'x'
	    	println("Old maximum radius (km) = ", arr[9])
    		www = inputf(" Enter new maximum radius (km) ")
		narr = (arr[1:8]..., www, arr[10:end]...)
		

	elseif stin == 'e'
    		println("Old nominal fractional accuracy for frequency = ", arr[10])
    		www = inputf(" Enter new nominal fractional accuracy ")
		narr = (arr[1:9]..., www, arr[11:end]...)

	elseif stin == 'd'
    		println("Previous number of dispersion curve points = ", arr[11])
    		iii = inputi(" Enter the new number of points on the dispersion curve ")
    		println("Old first azimuthal wavenumber, increment along curve (integers) = ", arr[12:13])
    		www = inputi(" Enter new first azimuthal wavenumber (integer) ")
    		if iii < 1.1
        		ar1 = arr[13]
    		else
        		ar1 = inputi(" Enter new azimuthal wavenumber increment (integer) ")
    		end
		

		narr = (arr[1:10]..., iii, www, ar1, arr[14:end]...)

	elseif stin == 'h'
    		println("Old number of depth input points = ", arr[14])
    		nd = inputi(" Enter number of (r, h) pairs  (> 0) ")
    		if nd < 1
        		println("Error!! nd must be >= 1")
			return
    		end
    		println("Old array of r locations for depth (km) = ", arr[15])	
    		ar1 = inputf(" Enter array of r locations where depth will be given (km) ")

    		println("Old array of depths (m) = ", arr[16])
    		ar2 = inputf(" Enter array of depths at locations given (m) ")
		
		iqq = findall(ar2 -> ar2 < 0, ar2)
		if length(iqq) > 0.5
			println(' ')
			println("Error! All depths must be > 0!")
			return 
		end

		if ar1[1] < 0.000001
			println(' ')
			println("	Error! The island must have a radius > 0!")
			return
		end
		if ar1[end] < arr[9]
			nd = nd + 1
			ar1 = [ar1 arr[9]]
			ar2 = [ar2 ar2[end]]
		end
		narr = (arr[1:13]..., nd, ar1, ar2, arr[17:end]...)


	elseif  stin == 'R'

    		println("Old number of R values = ", arr[17])
    		nR = inputi(" Enter number of (x, R) pairs (>=0) ")
  		if nR !=0
        	        println("Old R locations = ", arr[18])
        		ar1 = inputf(" Enter array of r locations where resistance coefficient is given (km) ")
            		println("Old R values = ", arr[19])
			ar2 = inputf(" Enter array of resistance coefficient (cm/sec) ")


			if length(ar1)!= nR
				println(' ')
				println("Error! Array length not equal to nR!")
				return
			end
			if length(ar2)!= nR
				println(' ')
				println("Error! Array length not equal to nR")
				return
			end
			radh = arr[15]
			if ar1[1] > radh[1]
				nR = nR + 1
				ar1 = [radh[1] ar1]
				ar2 = [ar2[1] ar2]
			end
			if ar1[end] < arr[9]
				nR = nR + 1
				ar1 = [ar1 arr[9]]
				ar2 = [ar2 ar2[end]]
			end

        	else
			nR = 2
			ar1 = [0 arr[9]]
			ar2 = [0 0]
		end
		narr = (arr[1:16]..., nR, ar1, ar2, arr[20:end]...)		
		
       
	elseif stin =='n'
		println("Old # of N^2 points = ", arr[20]) 
      		nn = inputi(" Enter new number of Nsquared values to read (>=1) ")
		println("Old depth increment for reading N^2 (m) = ", arr[21])
		ar1 = inputf(" Enter new depth increment for reading N^2 (m) ")
		println("Old exponential tail for need N^2 profile (km) = ", arr[22])
		ar2 = inputf(" Enter new exponential tail for N^2 (km) ")
		println("Old array of N^2 values (sec^-2) = ", arr[23])
		ar3 = inputf(" Enter new array of N^2 values (sec^-2) ")
		nn = length(ar3)
		if length(ar3) != nn
			println(' ')
			println("Error! Array length not equal to nn!")
			return
		end
		narr = (arr[1:19]..., nn, ar1, ar2, ar3, arr[24:end]...)
    
	elseif stin == 'v'

		println("Old azimuthal flow amplitude (cm/sec) = ", arr[24])	
 		va = inputf(" Enter mean azimuthal flow amplitude (cm/sec) ")
		if abs(va) > 0.0001
			println("Old radial distance of speed maximum (km) ", arr[25])
		        ar1 = inputf(" Enter new radial distance of velocity extreme (km) " )   
			println("Old vertical location of velocity extreme (m) = ", arr[26])
			ar2 = inputf(" Enter new vertical location of velocity extreme (m) ")
			println("Old downward and upward scales for mean velocity (m) ", arr[27], ",   ", arr[28])
			ar3 = inputf(" Enter new downward scale for mean flow (m) ")
			ar4 = inputf(" Enter new upward scale for mean flow ")
			println("Old offshore and onshore scales for mean flow (km) ", arr[29], ",  ", arr[30])
			ar5 = inputf(" Enter new offshore scale for mean velocity (km) ")
			ar6 = inputf(" Enter new onshore scale for mean velocity (km) ")
			if arr[31] < 0.5
				println("Old: Stratification is undisturbed at the center")
			else
				println("Old: Stratification is undisturbed far offshore")
			end
			istrat = inputi(" Enter 0 for stratification undisturbed at the center, 1 for undisturbed offshore ")
		else
			ar1 = 30.
			ar2 = 0.
			ar3 = 1000.
			ar4 = 1000.
			ar5 = 1000.
			ar6 = 1000.
			istrat = 1.
		end
		narr = (arr[1:23]..., va, ar1, ar2, ar3, ar4, ar5, ar6, istrat, arr[32]...)
    
	elseif stin =='p'
    	
    		if arr[end] > 0.5
		        println("Was set to pause for graphics")
		else
        		println("Was set not to pause for graphics")
    		end
    		ipause = inputi(" Enter 1 to pause for graphics or 0 to skip pauses ")

    		narr = (arr[1:31]..., ipause)
    
	else
    		println("Value entered does not correspond to any of the above categories")
		return
	end
	return narr

end									# end of function








"""
    inputi(quer)
    
    Replace the Matlab "input" function for integer numbers/arrays
  
  	quar = a string that asks for numbers/arrays

  It asks for numbers. You can either:
	type in a single number
	type in an array, e.g.:  [1 2 3]
     or (if you have the numbers in an array or vector x)
	type in string(x)
	
  
  Returns an integer number or an array of integers
  
  K. Brink  10/21/2025
"""

function inputi(quer)
	local respp, iar, icount, contt


    	print(quer)
    										#respp is a string
    	respp = readline()

	if length(respp) >= 6
		if respp[1:6] == "string"					# This allows you to enter an array by name
			ccc = Meta.parse(respp)
			rrr = eval(ccc)
			respp = rrr
		end
	end
    
    	while respp[1] == ' '            					# Remove any initial blanks
        	respp = respp[2:end]
    	end
    	while respp[end] == ' '	 						# remove any blanks at the end
		respp = respp[1:(end-1)]
	end
	
	
	if occursin(",",respp)							# remove commas and replace with blanks
		
		respp = replace(respp,',' => ' ')
	end
	
    
    	if respp[1] != '['							# single number, no brackets
        	iar = parse(Int,respp)
        	return iar
    	else									# input is a vector, not a scalar
		if respp[end] != ']'
			println("Error: missing bracket!")
			return NaN
		else
			respp = respp[2:(end-1)]    		 		#remove brackets
			while respp[end] ==' '
				respp = respp[1:(end-1)]
			end

			icount = 1;
			contt = 1;
			while contt > 0.5
				if occursin(' ',respp)
					ii = findfirst(' ',respp)
					tem = respp[1:(ii-1)];
					respp = respp[(ii+1):end]		# Shorten string
					while respp[1] == ' '
						respp = respp[2:end]
				 	end
				 	iii = parse(Int64,tem)
				 	if icount == 1
						iar = iii
				 	else
						iar = [iar iii]
				 	end
					icount = icount + 1
				else						# You are at the last number
					contt = 0				# signal that you are done
					iii = parse(Int64,respp)
					if icount > 1.5
						iar = [iar iii]
					else
						iar = iii
					end
				end
			end					
		
		return iar
		end						
	end							
	
end										# function end




"""
    inputf(quer)
    
  Replace the Matlab "input" function for floating point numbers/arrays
  
  	quer = a string that asks for numbers/arrays, e.g., 'Array of x (km)? "

  It asks for numbers. You can either:
	type in a single number
	type in an array, e.g.:  [1.1 2.6 3.2]
     or (if you have the numbers in an array or vector x)
	type in string(x)
  
  Returns a number or an array of floating point numbers
  
  K. Brink  10/21/2025
"""

function inputf(quer)
	local respp, flar, icount, contt
    	print(quer)
    
   	respp = readline()						# this is a string

	if length(respp) >= 6
		if respp[1:6] == "string"				# this allows you to enter an array by name
			ccc = Meta.parse(respp)
			rrr = eval(ccc)
			respp = rrr
		end
	end
    
    	while respp[1] == ' '            				# Remove any initial blanks
        	respp = respp[2:end]
    	end
	while respp[end] == ' '						# remove any blanks at the end
		respp = respp[1:(end-1)]
	end
	
	if occursin(",",respp)						# remove commas and replace with blanks
		
		respp = replace(respp,',' => ' ')
	end
	
		
	
    
    	if respp[1] != '['						# single number, no brackets
        	flar = parse(Float64,respp)
        	return flar
    	else								# input is a vector, not a scalar
		if respp[end] != ']'
			println("Error: missing bracket!")
			return NaN
		else
			respp = respp[2:(end-1)]   	 	 	#remove brackets
			while respp[end] ==' '
				respp = respp[1:(end-1)]
			end

			icount = 1
			contt = 1
			while contt > 0.5
				if occursin(' ',respp)
					ii = findfirst(' ',respp)
					tem = respp[1:(ii-1)]
					respp = respp[(ii+1):end]		# Shorten string
					while respp[1] == ' '
						respp = respp[2:end]
				    	end
					fii = parse(Float64,tem)
					if icount == 1
						flar = fii
					else
						flar = [flar fii]
					end
					icount = icount + 1
				else						# You are at the last number
					contt = 0				# signal that you are done
					fii = parse(Float64,respp)
					if icount > 1.5
						flar = [flar fii]
					else
						flar = fii
					end
				end
			end
			return flar
		end
	end
end										# function end


