FUNCTION MYFUNCTION, X, Parameters

voigt_profile_1 =  Parameters[0] * cnb_voigt(x-Parameters[1], Parameters[2], Parameters[3]) 

voigt_profile_2 =  Parameters[4] * cnb_voigt(x-Parameters[5], Parameters[6], Parameters[7]) 

line = Parameters[8] + x*Parameters[9]

RETURN, voigt_profile_1 + voigt_profile_2 + line 

END

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function month2number, month

month=STRLOWCASE(month)

CASE month OF
   'enero': return, 1
   'febrero': return, 2
   'marzo': return, 3
   'abril': return, 4
   'mayo': return, 5
   'junio': return, 6
   'julio': return, 7
   'agosto': return, 8
   'septiembre': return, 9
   'octubre': return, 10
   'noviembre': return, 11
   'diciembre': return, 12
   ELSE: stop
ENDCASE

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO fit_2_voigt

readcol, "names.list", name_spec, f="A"

nspec=n_elements(name_spec)

cgps_open, "testps_voigt.ps", scale_factor=1

;File where we are going to store the results:
text_file_name = "ew_voigt.txt"

;We check if the file already exist
if file_test(text_file_name) then begin

   repeat_counter = 0
   repeat begin

      answer = ""  
      READ, answer, PROMPT='The file "' + text_file_name + '" already exist. Do you want to overwrite it? (y/n)' 
      answer = STRLOWCASE(answer)
   
      ;We will ask ~3 times
      repeat_counter = repeat_counter + 1
      if repeat_counter gt 3 then return
      if answer eq "n" then return

   ;We force the answer to be "y" to continue
   endrep until answer eq "y"
   
endif

;Header of the file
openw, 1, text_file_name, width=300
printf, 1, "# EW_1"
printf, 1, "# EW_2"
printf, 1, "# Delta_EW_1"
printf, 1, "# Delta_EW_2"
printf, 1, "# JD"
printf, 1, "# FWHM_1"
printf, 1, "# FWHM_2"
printf, 1, "# Centroid_1"
printf, 1, "# Centroid_2"
printf, 1, "# Ruido_por_pixel"
printf, 1, "# Star_name"
printf, 1, "# plate"
close, 1

for ispec=0,nspec-1 do begin

   ; In case we want only one trial
   if n_elements(uno) ne 0 then ispec=uno
   
   ; Antes de graficar, hayaremos la fecha de los espectros
   readcol, "../../../../../Log_placas.txt", log_name, log_plate, log_day, log_month, log_year, log_prism, format="A,A,I,A,I,I", skip=1 ; <------------------
   star_name = ( strsplit(name_spec[ispec],".",/extract) )[0]
   plate = ( strsplit(star_name,"_",/extract) )[1]
   star_name = ( strsplit(star_name,"_",/extract) )[0]

   iwhere = where(log_name eq star_name and log_plate eq plate, nwhere)
   if nwhere ne 1 then stop

   title = star_name + " - plate " + plate + " - He $\lambda$ 4471 $\angstrom$ & Mg $\lambda$ 4481 $\angstrom$"; <---------------------
   date = string(log_day[iwhere]) + " " + log_month[iwhere] + " " + strtrim(string(log_year[iwhere]),2)

   ;Leemos el archivo fits.
   myreadfits, "../../../../../"+name_spec[ispec], lambda, flux_norm

   ; Esta es la region donde calcularemos el "ruido":  <---------------------
   ;;i=where(lambda ge 4600 and lambda le 4800)
   ;;lambda_ruido=lambda[i]
   ;;flux_ruido=flux_norm[i]

   ; Cortamos el espectro en torno a la linea.  <---------------------
   lambda_windows=[4240,4440]; <---------------------
   i=where(lambda ge lambda_windows[0] and lambda le lambda_windows[1])
   lambda=lambda[i]
   flux_norm=flux_norm[i]

   ;Apuesta inicial para MPFIT (obtenida APLICANDO MPFITPEAK sin restricciones) <---------------------
   start_params = [-1.541091, 4341.00, 4.21132, 4.21132, $
                   -1.541091, 4387.00, 4.21132, 4.21132, $
                   1.0, 0.0]

   ;Sin pesos
   err=replicate(1.,n_elements(lambda))
   weights=replicate(1.,n_elements(lambda))
   ;Lugares que queremos omitir en todos los espectros:
   ;;;weights[where(lambda ge 3927 and lambda le 3940)]=0 ;  <---------------------

   ;Notemos que los valores de "start_params" se usaran en lugar de los de .value <---------------------
   parinfo = replicate({value:0.D, fixed:0, limited:[0,0], limits:[0.,0], parname:"",tied:""}, 10)

   ;Notemos que los valores de "start_params" se usaran en lugar de los de .value <---------------------
   parinfo[*].value = start_params

   ;Names
   parinfo[*].PARNAME = ["Peak Value 1","Peak Centroid 1","Gaussian Sigma 1","Lorentz Sigma 1",$
                         "Peak Value 2","Peak Centroid 2","Gaussian Sigma 2","Lorentz Sigma 2",$
                         "Constant", "slope"]

   ;Slope and base
   parinfo[8].fixed = 1
   parinfo[9].fixed = 1
   
   ;position
   parinfo[1].limited[*] = 1
   parinfo[1].limits = [4330.0,4350.0]
   parinfo[5].limited[*] = 1
   parinfo[5].limits = [4383.0,4391.0]
   
   ;peak
   parinfo[0].limited[1] = 1
   parinfo[0].limits = 0.
   parinfo[4].limited[1] = 1
   parinfo[4].limits = 0.
   
   ;fwhm
   parinfo[2].limited[*] = 1
   parinfo[2].limits[*] = [0.,10.]
   parinfo[3].limited[*] = 1
   parinfo[3].limits[*] = [0.,10.]   
   parinfo[6].limited[*] = 1
   parinfo[6].limits[*] = [1.,10.]
   parinfo[7].limited[*] = 1
   parinfo[7].limits[*] = [1.,10.]
   

   ;Particular cases, lugares que queremos omitir:
   ;particular_case = 0
;;   if plate eq 1524 and star_name eq "HR4037" then begin
;;     weights[where(lambda ge 4390 and lambda le 4410)] = 0.
;;   endif



   ;fit  
   ;;;yfit = mpfitpeak(lambda, flux_norm, A, /lorentz, nterms=5, parinfo=parinfo, estimates=start_params, weights=weights)
   A = MPFITFUN('myfunction', lambda, flux_norm, err, start_params, yfit=yfit, parinfo=parinfo, STATUS=STATUS, weights=weights)

   line = A[8] + lambda*A[9]
   voigt_profile_1 =  A[0] * cnb_voigt(lambda-A[1], A[2], A[3]) 
   voigt_profile_2 =  A[4] * cnb_voigt(lambda-A[5], A[6], A[7]) 

   FWHM_1 = mean(minmax(voigt_profile_1 + line))
   FWHM_1 = minmax(lambda[where(voigt_profile_1 + line le FWHM_1)])
   FWHM_1 = FWHM_1[1] - FWHM_1[0]
   
   FWHM_2 = mean(minmax(voigt_profile_2 + line))
   FWHM_2 = minmax(lambda[where(voigt_profile_2 + line le FWHM_2)])
   FWHM_2 = FWHM_2[1] - FWHM_2[0]

   ; Graficamos
   cgplot, lambda, flux_norm, position=[0.231038,0.125,0.793962,0.9], xtickformat="(A1)", charsize=1.0, xr=lambda_windows, yr=[0,1.2];, asp=0.774 ; <---------------
   cgoplot, lambda, yfit, color="blue", thick=2   
   cgoplot, lambda, voigt_profile_1 + line, color="red", line=1   
   cgoplot, lambda, voigt_profile_2 + line, color="red", line=1   

   ;Graficamos tambien el continuo determinado en esta ventana
   cgoplot, lambda, A[8]+A[9]*lambda, line=2, color="blue"

   ;EWs. Notemos que ahora solo integramos al rededor de la linea de absorcion.
   ;if abs(A[0]) ge 0.001 then EW_range_1=[ A[1]-1.5*abs(FWHM_1), A[1]+1.5*abs(FWHM_1) ] else EW_range_1=[0,0]
   ;if abs(A[4]) ge 0.001 then EW_range_2=[ A[5]-1.5*abs(FWHM_2), A[5]+1.5*abs(FWHM_2) ] else EW_range_2=[0,0]
   if abs(A[0]) ge 0.001 then begin
      test_1 = where(voigt_profile_1 le -0.01, n_test_1)
      if n_test_1 gt 2 then EW_range_1 = minmax(lambda[where(voigt_profile_1 le -0.01)]) else EW_range_1=[0,0] 
   endif else EW_range_1=[0,0];EW_range_1=[ A[1]-1.5*abs(FWHM_1), A[1]+1.5*abs(FWHM_1) ]
   if abs(A[4]) ge 0.001 then begin
      test_2 = where(voigt_profile_2 le -0.01, n_test_2)
      if n_test_2 gt 2 then EW_range_2 = minmax(lambda[where(voigt_profile_2 le -0.01)]) else EW_range_2=[0,0];EW_range_2=[ A[5]-1.5*abs(FWHM_2), A[5]+1.5*abs(FWHM_2) ] 
   endif else EW_range_2=[0,0]   
   if EW_range_1[0] ne 0 then i_lambda_EW_1=where(lambda ge EW_range_1[0] and lambda le EW_range_1[1], ni_lambda_EW_1) else i_lambda_EW_1=-1
   if EW_range_2[0] ne 0 then i_lambda_EW_2=where(lambda ge EW_range_2[0] and lambda le EW_range_2[1], ni_lambda_EW_2) else i_lambda_EW_2=-1
   ;Notemos que vamos a usar el continuo calculado por el ajuste (parametros A[3] y A[4])
   ;EW_1 = INT_TABULATED( lambda[i_lambda_EW], yfit[i_lambda_EW]-1 + (A[3] + lambda[i_lambda_EW]*A[4]) )
   if i_lambda_EW_1[0] ne -1 then EW_1 = INT_TABULATED(lambda[i_lambda_EW_1], voigt_profile_1[i_lambda_EW_1]) else EW_1=0
   if i_lambda_EW_2[0] ne -1 then EW_2 = INT_TABULATED(lambda[i_lambda_EW_2], voigt_profile_2[i_lambda_EW_2]) else EW_2=0
   ;EW_1 = INT_TABULATED(lambda, yfit-1)

         
   ;Region donde calcularemos el EW y por ende el ruido. 
   if i_lambda_EW_1[0] ne -1 then cgoplot, EW_range_1, [1.10,1.10], color="blue", ps=-1;thick=2
   if i_lambda_EW_1[0] ne -1 then cgoplot, mean(EW_range_1), [1.10], color="green", ps=1;thick=2
   if i_lambda_EW_2[0] ne -1 then cgoplot, EW_range_2, [1.13,1.13], color="blue", ps=-1;thick=2
   if i_lambda_EW_2[0] ne -1 then cgoplot, mean(EW_range_2), [1.13], color="green", ps=1;thick=2

   
   ;Particular case, plots:
   if total(weights eq 0) gt 0 then begin
      i_lambda_particular_cases = where(abs(weights-shift(weights,1)) eq 1)
      vline, lambda[i_lambda_particular_cases], color="green"
      for k_i_lambda_particular_cases=0, n_elements(i_lambda_particular_cases) - 1, 2 do cgoplot, [ lambda[i_lambda_particular_cases[k_i_lambda_particular_cases]], lambda[i_lambda_particular_cases[k_i_lambda_particular_cases+1]] ], !y.crange, color="green"
   endif 

   ;Diferencia
   ;;cgplot, lambda, flux_norm-1-yfit, yr=[-0.1,0.1], position=[0.231038,0.0,0.793962,0.125], /noerase, xtickformat="(A1)", charsize=0.5
   !x.ticklen=0.1
   cgplot, lambda, flux_norm - yfit, yr=[-0.1,0.1], position=[0.231038,0.0,0.793962,0.125],/noerase, ytickformat="(A1)", xr=lambda_windows, charsize=1.0 ; <---------------
   !x.ticklen=0.0
   cgaxis, yaxis=0, yrange=[-0.1,0.1], charsize=1.0, color="blue", yTICKN=[' ', '-0.05', '0.00', '0.05', ' ']; ystyle=0
   cgaxis, yaxis=1, yrange=[-0.1,0.1], charsize=1.0, color="blue", yTICKN=[' ', ' ', ' ', ' ', ' ']; ystyle=0

   ;Particular case, plots:
   if total(weights eq 0) gt 0 then begin
      i_lambda_particular_cases = where(abs(weights-shift(weights,1)) eq 1)
      vline, lambda[i_lambda_particular_cases], color="green"
      for k_i_lambda_particular_cases=0, n_elements(i_lambda_particular_cases) - 1, 2 do cgoplot, [ lambda[i_lambda_particular_cases[k_i_lambda_particular_cases]], lambda[i_lambda_particular_cases[k_i_lambda_particular_cases+1]] ], !y.crange, color="green"
   endif 

   ;"Ruido"
   ;Plot de la region donde tomamos / calculamos el ruido
   ;cgplot, lambda_ruido, flux_ruido-1, yr=[-0.1,0.1], position=[0.0,0.9,1.0,1.0], /noerase, charsize=0.5, xtickformat="(A1)";, xstyle=4
   ;0.231038,0.9,0.793962,1.0
   ;;cgaxis, xaxis=1, xrange=[4600,4800], xstyle=1, charsize=0.5
   ;;;;;ruido_pix = stddev(flux_ruido) ;;; <------------
   ;;;;;ruido_pix = stddev(flux_norm[i_lambda_EW]-yfit[i_lambda_EW])
   ;;;;;ruido_pix = stddev(flux_norm-yfit)
   ;Tomamos la desviacion estandar en toda la ventana. Omitimos los casos particulares
   ruido_pix = stddev(flux_norm[where(weights eq 1)]-yfit[where(weights eq 1)])
   ;;cgtext, 1.00, 1.05,  "Ruido por pixel = "+string(ruido_pix,format="(G0)"), /normal, charsize=1.5, color="blue", ALIGNMENT=1.0
   cgtext, 0.5, 1.02,  "Ruido por pixel = "+string(ruido_pix,format="(G0)")+". Calculado como la desviacion estandar del residuo.", /normal, charsize=1.5, color="blue", ALIGNMENT=0.5
   
   ;EWs
   ;;;ruido_1 = sqrt(4*A[2])*ruido_pix
   if i_lambda_EW_1[0] ne -1 then ruido_1 = sqrt(ni_lambda_EW_1)*ruido_pix else ruido_1 = 0
   if i_lambda_EW_2[0] ne -1 then ruido_2 = sqrt(ni_lambda_EW_2)*ruido_pix else ruido_2 = 0

   ;Text
   cgtext, 0.25,0.35, "EW_1 = "+string(EW_1,format="(G0)"), /normal, charsize=1, color="red"
   cgtext, 0.25,0.30, "$\Delta$ = "+string(ruido_1,format="(G0)")+$
     "   %$\Delta$ = "+string(abs(ruido_1/EW_1*100),format="(G0)"), /normal, charsize=1, color="red"
   cgtext, 0.25,0.20, "EW_2 = "+string(EW_2,format="(G0)"), /normal, charsize=1, color="red"
   cgtext, 0.25,0.15, "$\Delta$ = "+string(ruido_2,format="(G0)")+$
     "   %$\Delta$ = "+string(abs(ruido_2/EW_2*100),format="(G0)"), /normal, charsize=1, color="red"

   ;Titulo
   cgtext, 0.50, 1.1,  title, /normal, charsize=2, ALIGNMENT=0.5

   ;Julian day
   JD = julday(month2number(log_month[iwhere]), log_day[iwhere], log_year[iwhere])

   openw, 1, text_file_name, /append, width=300
   printf, 1, EW_1, EW_2, ruido_1, ruido_2, JD[0], FWHM_1, FWHM_2, A[1], A[5], ruido_pix,  string(star_name, format="(A9)"), string(plate, format="(A9)")
   close, 1

   if n_elements(uno) ne 0 then break

endfor

;Parametros en el  PDF
cgplot, [1], [1], /nodata, ysty=4, xsty=4

;Paramerto [0]
cgtext, 0.00, 1.00, parinfo[0].parname, /normal, charsize=1, color="red"
cgtext, 0.00, 0.95, "starting value = "+string(start_params[0],format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.90, "fixed = "+string(parinfo[0].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.85, "limited = "+string(parinfo[0].limited[0],format="(G0)")+" "+$
                                 string(parinfo[0].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.80, "limits = "+string(parinfo[0].limits[0],format="(G0)")+" "+$
                                string(parinfo[0].limits[1],format="(G0)"), /normal, charsize=1
;Paramerto [1]
cgtext, 0.25, 1.00, parinfo[1].parname, /normal, charsize=1, color="red"
cgtext, 0.25, 0.95, "starting value = "+string(start_params[1],format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.90, "fixed = "+string(parinfo[1].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.85, "limited = "+string(parinfo[1].limited[0],format="(G0)")+" "+$
                                 string(parinfo[1].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.80, "limits = "+string(parinfo[1].limits[0],format="(G0)")+" "+$
                                string(parinfo[1].limits[1],format="(G0)"), /normal, charsize=1
;Paramerto [2]
cgtext, 0.50, 1.00, parinfo[2].parname, /normal, charsize=1, color="red"
cgtext, 0.50, 0.95, "starting value = "+string(start_params[2],format="(G0)"), /normal, charsize=1
cgtext, 0.50, 0.90, "fixed = "+string(parinfo[2].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.50, 0.85, "limited = "+string(parinfo[2].limited[0],format="(G0)")+" "+$
                                 string(parinfo[2].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.50, 0.80, "limits = "+string(parinfo[2].limits[0],format="(G0)")+" "+$
                                string(parinfo[2].limits[1],format="(G0)"), /normal, charsize=1

;Paramerto [3]
cgtext, 0.75, 1.00, parinfo[3].parname, /normal, charsize=1, color="red"
cgtext, 0.75, 0.95, "starting value = "+string(start_params[3],format="(G0)"), /normal, charsize=1
cgtext, 0.75, 0.90, "fixed = "+string(parinfo[3].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.75, 0.85, "limited = "+string(parinfo[3].limited[0],format="(G0)")+" "+$
                                 string(parinfo[3].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.75, 0.80, "limits = "+string(parinfo[3].limits[0],format="(G0)")+" "+$
                                string(parinfo[3].limits[1],format="(G0)"), /normal, charsize=1
                                
                                
                                
                                
                                

;Paramerto [4]
cgtext, 0.00, 0.75, parinfo[4].parname, /normal, charsize=1, color="red"
cgtext, 0.00, 0.70, "starting value = "+string(start_params[4],format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.65, "fixed = "+string(parinfo[4].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.60, "limited = "+string(parinfo[4].limited[0],format="(G0)")+" "+$
                                 string(parinfo[4].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.55, "limits = "+string(parinfo[4].limits[0],format="(G0)")+" "+$
                                string(parinfo[4].limits[1],format="(G0)"), /normal, charsize=1
;Paramerto [5]
cgtext, 0.25, 0.75, parinfo[5].parname, /normal, charsize=1, color="red"
cgtext, 0.25, 0.70, "starting value = "+string(start_params[5],format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.65, "fixed = "+string(parinfo[5].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.60, "limited = "+string(parinfo[5].limited[0],format="(G0)")+" "+$
                                 string(parinfo[5].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.55, "limits = "+string(parinfo[5].limits[0],format="(G0)")+" "+$
                                string(parinfo[5].limits[1],format="(G0)"), /normal, charsize=1
                                
;Paramerto [6]
cgtext, 0.50, 0.75, parinfo[6].parname, /normal, charsize=1, color="red"
cgtext, 0.50, 0.70, "starting value = "+string(start_params[6],format="(G0)"), /normal, charsize=1
cgtext, 0.50, 0.65, "fixed = "+string(parinfo[6].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.50, 0.60, "limited = "+string(parinfo[6].limited[0],format="(G0)")+" "+$
                                 string(parinfo[6].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.50, 0.55, "limits = "+string(parinfo[6].limits[0],format="(G0)")+" "+$
                                       string(parinfo[6].limits[1],format="(G0)"), /normal, charsize=1                         

;Paramerto [7]
cgtext, 0.75, 0.75, parinfo[7].parname, /normal, charsize=1, color="red"
cgtext, 0.75, 0.70, "starting value = "+string(start_params[7],format="(G0)"), /normal, charsize=1
cgtext, 0.75, 0.65, "fixed = "+string(parinfo[7].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.75, 0.60, "limited = "+string(parinfo[7].limited[0],format="(G0)")+" "+$
                                 string(parinfo[7].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.75, 0.55, "limits = "+string(parinfo[7].limits[0],format="(G0)")+" "+$
                                string(parinfo[7].limits[1],format="(G0)"), /normal, charsize=1
                                
                                
                                
                                
                                
                                
;Paramerto [8]
cgtext, 0.00, 0.50, parinfo[8].parname, /normal, charsize=1, color="red"
cgtext, 0.00, 0.45, "starting value = "+string(start_params[8],format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.40, "fixed = "+string(parinfo[8].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.35, "limited = "+string(parinfo[8].limited[0],format="(G0)")+" "+$
                                 string(parinfo[8].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.00, 0.30, "limits = "+string(parinfo[8].limits[0],format="(G0)")+" "+$
                                string(parinfo[8].limits[1],format="(G0)"), /normal, charsize=1
;Paramerto [9]
cgtext, 0.25, 0.50, parinfo[9].parname, /normal, charsize=1, color="red"
cgtext, 0.25, 0.45, "starting value = "+string(start_params[9],format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.40, "fixed = "+string(parinfo[9].fixed,format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.35, "limited = "+string(parinfo[9].limited[0],format="(G0)")+" "+$
                                 string(parinfo[9].limited[1],format="(G0)"), /normal, charsize=1
cgtext, 0.25, 0.30, "limits = "+string(parinfo[9].limits[0],format="(G0)")+" "+$
                                string(parinfo[9].limits[1],format="(G0)"), /normal, charsize=1


cgps_close, /pdf

return

END
