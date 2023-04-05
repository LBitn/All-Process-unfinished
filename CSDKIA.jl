module CSDKIA

# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
GaussianKernel = [
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
    0.00002292	0.00078634	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
    0.00038771	0.01330373	0.11098164	0.22508352	0.11098164	0.01330373	0.00038771
    0.00019117	0.00655965	0.05472157	0.11098164	0.05472157	0.00655965	0.00019117
    0.00002292	0.00078633	0.00655965	0.01330373	0.00655965	0.00078633	0.00002292
    0.00000067	0.00002292	0.00019117	0.00038771	0.00019117	0.00002292	0.00000067
];
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #
# El operador de Laplace-Lindenberg
Laplacian01 = [ [ 0 1 0 ]; [ 1 -4 1 ]; [ 0 1 0 ] ];
Laplacian02 = [ [ 0.5 0 0.5 ]; [ 0 -2 0 ]; [ 0.5 0 0.5 ] ];
LaplacianKernel = ( ( 1 - ( 1 / 3 ) ) * Laplacian01 ) + ( ( 1 / 3 ) * Laplacian02 );
# •·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·•·•·•·••·•·•·•·•·•·• #

using StatsBase
using PyPlot

export AntesQue
export DespuesQue
export ImagenTrayectoriasPosyNeg
export ImagenTraysconCSDA
export ImagenDeMosaico
export TiraOrillas
export vecindad8
export ComponentesSP
export ObtenComponentesyCM
export dist2D
export encuentraTrayectorias
export CSDAK
export CSDAI
export DiscreteLaplacian
export GaussianSmooth
export GaussSuavizarTemporalI
export GaussSuavizarTemporalK
export UnNormGaussI
export UnNormGaussK

function TiraOrillas( Puntos::Set )
    #Descarta lo que se sale de la malla de electrodos
    result = Set( [ ] )
    for p in Puntos
        if !( p[ 1 ] == 0 || p[ 2 ] == 0 || p[ 1 ] == 65 || p[ 2 ] == 65 )
            push!( result, p );
        end
    end
    return result
end

function vecindad8( punto::Array )
    # La ocho-vecindad de un punto en una malla cuadrada.
    j = punto[ 1 ];
    k = punto[ 2 ];
    result = Set{ Array{ Int64, 1 } }( );
    push!( result, [ j - 1, k - 1 ] );
    push!( result, [ j - 1, k ] );
    push!( result, [ j - 1, k + 1 ] );
    push!( result, [ j, k - 1 ] );
    push!( result, [ j, k + 1 ] );
    push!( result, [ j + 1, k - 1 ] );
    push!( result, [ j + 1, k ] );
    push!( result, [ j + 1, k + 1 ] );
    result = TiraOrillas( result )
    return result
end

function ComponentesSP( DatosSignados::Array )
    #Single pass method for Disjoint Components.
    lista = copy( DatosSignados );
    componentes = Set{ Any }( );
    while ( length( lista ) != 0 )
        x = pop!( lista ); #arranca el ULTIMO elemento de la lista
        listaprofundeza = Array{ Int64 }[ ];
        componentecurlab = Array{ Int64 }[ ];
        push!( listaprofundeza, x ); #Pone elementos al FINAL de la lista
        push!( componentecurlab, x );
        profundidad = 0;
        while ( ( length( listaprofundeza ) != 0 ) && profundidad < 1000 )
            y = pop!( listaprofundeza );
            for v in vecindad8( y )
                if in( v, lista ) # A: Si v está en la lista
                    deleteat!( lista, indexin( Any[ v ], lista ) );
                    push!( listaprofundeza, v );
                    profundidad += 1;
                    push!( componentecurlab, v );
                end
            end
        end
        push!( componentes, componentecurlab );
    end
    return componentes
end

function ObtenComponentesyCM( Datos::Array, tini = 1, tfini = tmax, epsilon = 1.0 )
    ( alto, ancho, lu ) = size( Datos );
    tamano = 3;
    CMPositivo=Dict{Int, Array}()
    CMNegativo=Dict{Int, Array}()
    for t=tini:tfini
        ActividadNegativa=Array{Int16}[]
        ActividadPositiva=Array{Int16}[]
        SpikeCountPositivo=zeros(alto,ancho)
        SpikeCountNegativo=zeros(alto,ancho)
        for j=1:alto,k=1:ancho
            if(Datos[j,k,t]<-epsilon)
                push!(ActividadNegativa, [j, k])
                SpikeCountNegativo[j,k]+=1
            elseif(Datos[j,k,t]>epsilon)
                push!(ActividadPositiva, [j, k])
                SpikeCountPositivo[j,k]+=1
            end
        end
        componentesneg=ComponentesSP(ActividadNegativa)
        centrosdemasaneg=[[0 0 0];]
        for p in componentesneg
            mu=length(p)
            if mu>tamano
                masa=0.00
                x=0.00
                y=0.00
                for q in p
                    j=q[1]
                    k=q[2]
                    masalocal=Datos[j,k,t]
                    masa+=masalocal
                    x+=k*masalocal
                    y+=j*masalocal
                end
                x/=masa                    # A: Ahora x es igual a x/masa
                y/=masa
                A=[x y masa]               # A: A son los datos del centro de masa del conjunto actualvcat(centrosdemasaneg, A)
                centrosdemasaneg=vcat(centrosdemasaneg, A)
            end
        end
        centrosdemasaneg=centrosdemasaneg[2:end,:]
        CMNegativo[t]=centrosdemasaneg
        componentespos=ComponentesSP(ActividadPositiva)
        centrosdemasapos=[[0 0 0];]
        for p in componentespos
            mu=length(p)
            if mu>tamano
                masa=0.00
                x=0.00
                y=0.00
                for q in p
                    j=q[1]
                    k=q[2]
                    masalocal=Datos[j,k,t]
                    masa+=masalocal
                    x+=k*masalocal
                    y+=j*masalocal
                end
                x/=masa
                y/=masa
                A=[x y masa]
                centrosdemasapos=vcat(centrosdemasapos, A)
            end
        end
        centrosdemasapos=centrosdemasapos[2:end,:]
        CMPositivo[t]=centrosdemasapos
    end
    return (CMPositivo, CMNegativo)
end
function dist2D(x,y)
    result=sqrt((x[1]-y[1])^2+(x[2]-y[2])^2)
    return result
end
function encuentraTrayectorias( Datos, toleradist, mincadena=20, mingordo=2.0, desde=1,hasta=20, )
    tau=1
    t=1
    j=1
    Catenario=Dict{Integer, Array{Any}}()
    Cadena=[0 0 0 0]
    tnum=1
    CopiaMegaArray=deepcopy(Datos);
    NumFrames=length(Datos)
    FakeNumFrames=NumFrames
    while t <= FakeNumFrames-1
        tau=t
        @label arrrrh
            if(CopiaMegaArray[tau]==[])

                jmax,nada=0,0
            else
         jmax,nada= size(CopiaMegaArray[tau])
            end
        while j <=jmax && tau<FakeNumFrames
                if abs(CopiaMegaArray[tau][j,3]) > mingordo
                Eslabon=[transpose(CopiaMegaArray[tau][j,:]) tau]
                Cadena=vcat(Cadena, Eslabon)
                mindist=2
                kasterisco=1
                    if CopiaMegaArray[tau+1]==[]
                        kmax,nada=0,0
                    else
                    kmax, nada= size(CopiaMegaArray[tau+1])
                    end
                    huboalgo=false
                for k=1:kmax
                    EslabonTentativo=CopiaMegaArray[tau+1][k,:]
                        if abs(EslabonTentativo[3])>mingordo
                        dist=dist2D(Eslabon,EslabonTentativo)
                        if dist<mindist
                            mindist=dist
                            kasterisco=k
                            huboalgo=true
                        end
                    end
                end
                if huboalgo && mindist<toleradist
                    CopiaMegaArray[tau][j,3]=0.0000
                    if tau+1<FakeNumFrames
                        tau+=1
                        j=kasterisco
                        @goto arrrrh
                    else
                        Eslabon=[transpose(CopiaMegaArray[tau+1][kasterisco,:]) tau+1]
                        Cadena=vcat(Cadena, Eslabon)
                        j+=1
                        tau=t
                        if size(Cadena)[1]>mincadena
                            Catenario[tnum]=Cadena[2:end,:]
                            tnum+=1
                        end
                        Cadena=[0 0 0 0]
                        @goto arrrrh
                    end
                else
                    if size(Cadena)[1]>mincadena
                            Catenario[tnum]=Cadena[2:end,:]
                            tnum+=1
                    end
                    Cadena=[0 0 0 0]
                    j+=1
                    tau=t
                    @goto arrrrh
                end
            end
            j+=1
            tau=t
        end
        @label urrr
        j=1
        t+=1
        tau=t
        Cadena=[0 0 0 0]
    end
    return Catenario
end
function ImagenTrayectoriasPosyNeg( CatenarioPositivo , CatenarioNegativo )
    figure(figsize=(9,8))
    #axis("equal")
    xlim(0.0,65.0)
    ylim(0.0,65.0)
    tolerancia=0
    minlong=0
    maxlong=20000
    for p in values(CatenarioPositivo)
        gordura=abs(p[:3])
        longus,gordus=size(p)
        if (mean(gordura)>tolerancia) && (longus>minlong) && longus < maxlong
            xxpos=p[:,1]
            yypos=p[:,2]
            tiempos=p[:,4]/7022
            plot(xxpos, yypos, marker="o", markersize=0.25, color="r", lw=0.25, zorder=15)
            colores=scatter(xxpos,yypos, s=gordura*0.4, edgecolors="none",
            c=tiempos, cmap="autumn", vmin=0.0, vmax=0.72 )
            principios=scatter(xxpos[1],yypos[1], s=gordura*0.5,
            edgecolors="black",
            color="g", marker="s", label="Inicio" ,
            zorder=1)
            finales=scatter(xxpos[end],yypos[end], s=gordura*0.5, edgecolors="black",
            color="gold", marker="D", label="Final", zorder=2 )
        end
    end
    savefig( "TrayectoriasP.svg" , dpi = 92 )
    figure(figsize=(9,8))
    #axis("equal")
    xlim(0.0,65.0)
    ylim(0.0,65.0)
    for p in values(CatenarioNegativo)
        gordura=abs(p[:3])
        longus,gordus=size(p)
        if (mean(gordura)>tolerancia) && (longus>minlong) && longus < maxlong
            xxpos=p[:,1]
            yypos=p[:,2]
            tiempos=p[:,4]/7022
            plot(xxpos, yypos, marker="o", markersize=0.25, color="r", lw=0.25, zorder=15)
            colores=scatter(xxpos,yypos, s=gordura*0.4, edgecolors="none",
            c=tiempos, cmap="autumn", vmin=0.0, vmax=0.72 )
            principios=scatter(xxpos[1],yypos[1], s=gordura*0.5,
            edgecolors="black",
            color="g", marker="s", label="Inicio" ,
            zorder=1)
            finales=scatter(xxpos[end],yypos[end], s=gordura*0.5, edgecolors="black",
            color="gold", marker="D", label="Final", zorder=2 )
        end
    end
    savefig( "TrayectoriasN.svg" , dpi = 92 )
end
function ImagenTraysconCSDA( csda , CatenarioPositivo , CatenarioNegativo , frecuencia )
    escribevelo=false
    cuadro=1000
    retraso=0
    ImagenCSD=csda[:,:,cuadro];
    figure(figsize=(6,6))
    milisec=round((cuadro-retraso)/frecuencia,digits=1)
    title("t= $milisec ms")
    tick_params(labelbottom="on", labelleft="on", direction="out")
    PyPlot.xticks(fontsize=10)
    PyPlot.yticks(fontsize=10)
    PyPlot.xlim(0,65)
    PyPlot.ylim(0,65)
    limcsd=75
    bolitasrojas=0
    bolitasazules=0
    tight_layout()
    guacafondo=PyPlot.imshow(ImagenCSD, cmap="bwr", interpolation="bicubic",
    vmin=-limcsd, vmax=limcsd, extent=[1,65,1,65], origin="lower")
    for (k,p) in CatenarioPositivo
        paux=p
        paux=p
        cucho,fleto=size(paux)
        longus,falsus=size(p)
        if cucho>0
            #Si son chiquitos no nos interesan.
            xxpos=paux[:,1]
            yypos=paux[:,2]
            gordis=abs.(map(Float32, paux[:,3]))
            #println( gordis , length( gordis ))
            tiempos=round.(paux[:,4]/7.022,digits=1)
            #println( 0.06.* gordis , length( 0.06.* gordis ))
            PyPlot.plot(xxpos, yypos, marker="o", markersize=1, color="maroon", lw=1,zorder=1900)
            PyPlot.scatter(xxpos[end],yypos[end],marker="o", s=0.06* gordis[end] , color="red",alpha=0.3)    #
            inicios=PyPlot.scatter(xxpos[1],yypos[1], s=20, edgecolors="black",
            linewidth=1,
            facecolor="maroon", marker="o",zorder=1999 )
            finales=PyPlot.scatter(xxpos[end],yypos[end], s=15, edgecolors="orange",
            linewidth=1,facecolor="crimson", marker="^", label="Final",zorder=2000 )
            if escribevelo && (length(tiempos)>1)
                dist=norm([xxpos[end]-xxpos[end-1],yypos[end]-yypos[end-1]],2)
                vel=round(dist*frecuencia,1)
                annotate(vel, (xxpos[end], yypos[end]), fontsize=10, color="maroon", zorder=1950)
            end
        end #(cierra sobre cucho)
    end #cierra sobre cmpred
    for (k,p) in CatenarioNegativo
        paux=p
        paux=p
        cucho,fleto=size(paux)
        longus,falsus=size(p)
        if  cucho>0
            xxpos=paux[:,1]
            yypos=paux[:,2]
            gordis=abs.(map(Float32, paux[:,3]))
            tiempos= round.((paux[:,4].-retraso)/7.022,digits=1)
            PyPlot.plot(xxpos, yypos, marker="o", markersize=1, color="cyan", lw=1,zorder=1998)
            PyPlot.scatter(xxpos[end],yypos[end],marker="o", s=0.06*gordis[end], color="blue",alpha=0.3)
            inicios=PyPlot.scatter(xxpos[1],yypos[1], s=15, edgecolors="black",
            linewidth=1,
            facecolor="midnightblue", alpha=0.7, marker="o", label="Final",zorder=1999 )
            finales=PyPlot.scatter(xxpos[end],yypos[end], s=15, edgecolors="cyan",
            linewidth=1,
            facecolor="midnightblue", marker="^", label="Final",zorder=2000 )
            if(escribevelo) && (length(tiempos)>1)
                dist=norm([xxpos[end]-xxpos[end-1],yypos[end]-yypos[end-1]],2)
                vel=round(dist*frecuencia,1)
                annotate(vel, (xxpos[end], yypos[end]), fontsize=10, color="cyan", zorder=1950)
            end
        end #sobre cucho
    end #sobre cmnegred
    savefig( "TrayectCSD.svg" , dpi = 90 )
end

function ImagenDeMosaico( csda , CatenarioPositivo , CatenarioNegativo, frecuencia )
    ion()
    farofa, lista=subplots(5,5, figsize=(25,20))
    limcsd=400
    inicio=1084
    fin=1500
    retraso=inicio
    paso=14
    epsi=9
    vepsi=(-epsi,epsi)
    cepsi=("teal","tomato")
    for j=1:25
        pu=ceil(Int,j/5)
        pa=mod(j,5)
        if pa==0; pa=5;end
        cuadro=inicio+j*paso
        milisec=round((cuadro-retraso)/frecuencia,digits = 1)
        subplots_adjust(bottom=0.06,left=0.1)
        lista[pu,pa][:imshow](csda[:,:,cuadro],cmap="seismic",vmin=-limcsd,vmax=limcsd,
        interpolation="nearest",extent=[1.5,64.5,63.5,1.5])
        lista[pu,pa][:contour](csda[:,:,cuadro], vepsi, colors=cepsi)
        lista[pu,pa][:set_title]("t = $milisec ms")
        lista[pu,pa][:tick_params](bottom="off", axis="both",which="both", labelbottom="off", labelleft="off")
        for (k,p) in CatenarioPositivo
            paux=AntesQue(p,cuadro+paso)
            paux=DespuesQue(paux,retraso-paso)
            cucho,fleto=size(paux)
            longus,falsus=size(p)
            if  (p[end,4] >= cuadro+paso) && cucho>0
                xxpos=paux[:,1]
                yypos=paux[:,2]
                gordis=abs.(map(Float32, paux[:,3]))
                lista[pu,pa][:set_xlim]([0,65])
                lista[pu,pa][:set_ylim]([0,65])
                lista[pu,pa][:plot](xxpos, yypos, marker="o", markersize=1, color="red", lw=1.5,zorder=1998)
                finales=lista[pu,pa][:scatter](xxpos[end],yypos[end], s=15, edgecolors="black",
                linewidth=1.5,
                facecolor="red", marker="o", label="Final" ,zorder=2000)

            end
        end
        for (k,p) in CatenarioNegativo
            paux=AntesQue(p,cuadro+paso)
            paux=DespuesQue(paux,retraso)
            cucho,fleto=size(paux)
            longus,falsus=size(p)
            if  p[end,4]>=cuadro+paso && cucho>0
                #Si son chiquitos no nos interesan.
                xxpos=paux[:,1]
                yypos=paux[:,2]
                gordis=abs.(map(Float32, paux[:,3]))
                lista[pu,pa][:set_xlim]([0,65])
                lista[pu,pa][:set_ylim]([0,65])
                lista[pu,pa][:plot](xxpos, yypos, marker="o", markersize=1, color="darkblue", lw=1.5,zorder=1998)
                finales=lista[pu,pa][:scatter](xxpos[end],yypos[end], s=15, edgecolors="black",
                linewidth=1.5,
                facecolor="darkblue", marker="o", label="Final" ,zorder=2000)
            end
        end
    end
    savefig( "Mosaico.svg" )
end
function AntesQue(Datos::Array, tiempo)
    Cadena=[0 0 0 0]
    for a in eachindex(Datos[:,4])
        if Datos[a,4]<tiempo
            Cadena=vcat(Cadena, transpose(Datos[a,:]))
        end
    end
    Cadena=Cadena[2:end,:]
    return Cadena
end
function DespuesQue(Datos::Array, tiempo)
    Cadena=[0 0 0 0]
    for a in eachindex(Datos[:,4])
        if Datos[a,4]>tiempo
            Cadena=vcat(Cadena, transpose(Datos[a,:]))
        end
    end
    Cadena=Cadena[2:end,:]
    return Cadena
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function UnNormGaussK( x, sigma )
    return exp( ( -x * x ) / ( 2 * sigma ) );
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function UnNormGaussI( x, sigma )
    cte = ( 1 / ( sigma * sqrt( 2 * pi ) ) );
    return cte * exp( -( ( x ^ 2 ) / ( 2 * ( sigma ^ 2 ) ) ) );
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function GaussSuavizarTemporalK( Datos, Sigma = 3 )
    # Un suavizado Gaussiano temporal.
    # Esto es escencialmente un filtro pasabajos.
    # Depende implicitamente de la frecuencia de muestreo.
    # sigma esta medido en pixeles, es la desviacion estandar de nuestro kernel.
    # El medioancho de nuestra ventana seran 3*sigma
    medioancho = ceil( Sigma * 3 );
    colchon = ones( medioancho );
    result = zeros( size( Datos ) );
    datoscolchon = vcat( colchon * Datos[ 1 ], Datos, colchon * Datos[ end ] );
    kernel = map( x -> UnNormGaussK( x, Sigma ), collect( -medioancho : medioancho ) );
    kernel = kernel / ( sum( kernel ) );
    #La convolucion asi normalizada preserva el valor RELATIVO entre los puntos de la funcion.
    for t = ( medioancho + 1 ) : ( length( Datos ) + medioancho )
        result[ t - medioancho ] = sum(
            datoscolchon[ ( t - medioancho ) : ( t + medioancho ) ] .* kernel );
    end
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function GaussSuavizarTemporalI( Datos, Sigma = 3 )
    # Un suavizado Gaussiano temporal.
    # Esto es escencialmente un filtro pasabajos.
    # Depende implicitamente de la frecuencia de muestreo.
    # sigma esta medido en pixeles, es la desviacion estandar de nuestro kernel.
    # El medioancho de nuestra ventana seran 3*sigma
    medioancho = ceil( Sigma * 3 );
    colchon = ones( medioancho );
    result = zeros( size( Datos ) );
    datoscolchon = vcat( colchon * Datos[ 1 ], Datos, colchon * Datos[ end ] );
    kernel = map( x -> UnNormGaussI( x, Sigma ), collect( -medioancho : medioancho ) );
    kernel = kernel / ( sum( kernel ) );
    #La convolucion asi normalizada preserva el valor RELATIVO entre los puntos de la funcion.
    for t = ( medioancho + 1 ) : ( length( Datos ) + medioancho )
        result[ t - medioancho ] = sum(
            datoscolchon[ ( t - medioancho ) : ( t + medioancho ) ] .* kernel );
    end
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function GaussianSmooth( temp )
    result = zeros( size( temp ) );
    ( mu, lu ) = size( temp );
    # Primero, hacemos el padding con copia de los datos para que no se suavice demasiado
    ## Okey, parece que los imbeciles de rioarriba cambiaron la sintaxis de
    # rebanadas de matriz. Ahora CUALQUIER rebanada de matriz es colvec.
    arriba = reshape( temp[ 1, : ], ( 1, lu ) );
    abajo = reshape( temp[ end, : ], ( 1, lu ) );
    arr3 = vcat( arriba, arriba, arriba );
    aba3 = vcat( abajo, abajo, abajo );
    temp = vcat( arr3, temp, aba3 );
    for j = 1 : 3
        temp = hcat( temp[ :, 1 ], temp, temp[ :, end ] );
    end
    for j = 4 : ( mu + 3 ), k = 4 : ( lu + 3 )
        # los indices van primero, "renglones", luego "columnas", etc
        aux = temp[ ( j - 3 ) : ( j + 3 ), ( k - 3 ) : ( k + 3 ) ];
        result[ ( j - 3 ), ( k - 3 ) ] = sum( GaussianKernel .* aux );
    end
    # Esta convolución no respeta norma L2
    # result = result*maximum(abs(Datos))/maximum(abs(result))
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function DiscreteLaplacian( temp )
    ( mu, lu ) = size( temp );
    izq = reshape( temp[ 1, : ], ( 1, lu ) );
    der = reshape( temp[ end, : ], ( 1, lu ) );
    # Primero, hacemos el padding con copia de los datos para que no se suavice demasiado
    temp = vcat( izq, temp, der );
    temp = hcat( temp[ :, 1 ], temp, temp[ :, end ] );
    largo, ancho = size( temp );
    aux = Array{ Float32 }( undef, 3, 3 );
    result = zeros( size( temp ) );
    # A: En esta parte lo que se hace es calcular el CSD aplicando el Kernel laplaciano
    # a cada celda más su 8-Vecindad y posteriormente suma todos los resultados como valor
    # de la celda
    for j = 2 : ( largo - 1 ), k = 2 : ( ancho - 1 )
        #los indices van primero, "renglones", luego "columnas", etc
        aux = temp[ ( j - 1 ) : ( j + 1 ), ( k - 1 ) : ( k + 1 ) ];
        result[ j, k ] = sum( LaplacianKernel .* aux );
    end
    # DO Crop the borders
    result = result[ 2 : ( end - 1 ) , 2 : ( end - 1 ) ]
    return result
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function CSDAK( Data )
    nChs, nFrs = size( Data );
    side = Int( sqrt( nChs ) );
    Data3D = reshape( Data, side, side, nFrs );
    ( μ, ν, ι )  = size( Data3D );
    # We apply a Temporal Gaussian smoothing ( this greatly affects the animations )
    Data3Plain = zeros( μ, ν, ι );
    for j = 1:μ, l = 1:ν
        channel = vec( Data3D[ j, l, : ] );
        Data3Plain[ j, l, : ] = GaussSuavizarTemporalK( channel );
    end
    Φ = zeros( μ, ν, ι );
    ∇ = zeros( μ, ν, ι );
    # We spatially smooth the LFP with a two-dimensional Gaussian filter.
    # Later we obtain the dCSD.
    for τ = 1:ι
        Φ[ :, :, τ ] = GaussianSmooth( Data3Plain[ :, :, τ ] );
        ∇[ :, :, τ ] = DiscreteLaplacian( Φ[ :, :, τ ] );
    end
    ∇ = -1 .* ∇;
    return ∇
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
function CSDAI( Data )
    nChs, nFrs = size( Data );
    side = Int( sqrt( nChs ) );
    Data3D = reshape( Data, side, side, nFrs );
    ( μ, ν, ι )  = size( Data3D );
    # We apply a Temporal Gaussian smoothing ( this greatly affects the animations )
    Data3Plain = zeros( μ, ν, ι );
    for j = 1:μ, l = 1:ν
        channel = vec( Data3D[ j, l, : ] );
        Data3Plain[ j, l, : ] = GaussSuavizarTemporalI( channel );
    end
    Φ = zeros( μ, ν, ι );
    ∇ = zeros( μ, ν, ι );
    # We spatially smooth the LFP with a two-dimensional Gaussian filter.
    # Later we obtain the dCSD.
    for τ = 1:ι
        Φ[ :, :, τ ] = GaussianSmooth( Data3Plain[ :, :, τ ] );
        ∇[ :, :, τ ] = DiscreteLaplacian( Φ[ :, :, τ ] );
    end
    ∇ = -1 .* ∇;
    return ∇
end
# ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ ~ ~~ #
end
