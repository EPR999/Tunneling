using LinearAlgebra
using GenericLinearAlgebra

const k = 3.0
const xf = 1.2
const xb = 1.0
const twopi = 2 * π
const e2 = k * xf * (exp(-8 * xb * (2 * xf - xb)) / (1 - exp(-32 * xf * xb)))
const stdrad = 0.1#初期面の半径
const dimension = 2

module KScatteringMap
        const k = 3.0
        const xf = 1.2
        const xb = 1.0
        const e2 = k * xf * (exp(-8 * xb * (2 * xf - xb)) / (1 - exp(-32 * xf * xb)))
        #正の時間方向の写像
	function dotV(x)
		x = BigFloat(x)
		return k * x * exp(-8 * x^2) - e2 * (exp(-8 * (x - xb)^2) - exp(-8 * (x + xb)^2))
	end

	function ddotV(x)
		x = BigFloat(x)
	    return k*(1 - 16*x^2) * exp(-8*x^2) - 16*e2*((x - xb) * exp(-8*(x - xb)^2) - (x + xb) * exp(-8*(x + xb)^2))
	end
	function U(qp)
		return [qp[1] + qp[2] - BigFloat(0.5) * dotV(qp[1]), qp[2] - BigFloat(0.5) * dotV(qp[1]) - BigFloat(0.5) * dotV(qp[1] + qp[2] - BigFloat(0.5) * dotV( qp[1] ) )]
    end
    #逆写像
    function Ui(qp)
            return [BigFloat(qp[1]) - BigFloat(qp[2]) - BigFloat(0.5) * BigFloat(dotV(qp[1])), BigFloat(qp[2]) + BigFloat(0.5) * BigFloat(dotV(qp[1])) + BigFloat(0.5) * BigFloat(dotV(qp[1])) - BigFloat(qp[2]) - BigFloat(0.5) * BigFloat(dotV(qp[1]))]
    end
end


function dotV(x)
	return k * x * exp(-8 * x^2) - e2 * (exp(-8 * (x - xb)^2) - exp(-8 * (x + xb)^2))
end

function ddotV(x)
	return k*(1 - 16*x^2) * exp(-8*x^2) - 16*e2*((x - xb) * exp(-8*(x - xb)^2) - (x + xb) * exp(-8*(x + xb)^2))
end

function strictJ(qp)
    DU = [ 1 - 0.5 * ddotV(qp[1]) , 1 ,- 0.5 * ddotV(qp[1]) - 0.5 * (1 - 0.5 * ddotV(qp[1])) * ddotV(qp[1] + qp[2] - 0.5 * dotV(qp[1])), 1 - BigFloat(0.5) * ddotV(qp[1] + qp[2]  - BigFloat(0.5) * dotV(qp[1]))  ]
	DU = reshape(DU,2,2)
	return DU
end

function strictJ2(func,qp,period)
	jacobian = strictJ(qp)
	for i = 1 : period-1
		qp = func(qp)
		jacobian = jacobian * strictJ(qp)
	end
	return jacobian
end

function NR(func,x,period)
	for i = 1:100
		x = x - LinearAlgebra.inv(Dg(func,x,period)) * g(func,x,period)
		if norm(g(func,x,period)) < 10 ^ (-18)
			println("break")
			return x
		elseif norm(x) > 3
				break
		end
	end
	return [100,100]
end

function randompp(period,func)#get from NR randomly if p-1 pp is insufficient to use for initial condition
	seedn = 300
	box = [-1.2,0.0]
    for i = 1 : seedn
		x = [BigFloat(rand()*2.6 - 1.3),BigFloat(rand()*2.6 - 1.3)]
        x = NR(func,x,period)
        #println(x)
		if norm(x) < 10
			box = hcat(box,[x[1],x[2]])
			box = del(box,[x[1],x[2]])
			println(box)
		end
	end
	println(box)
	#box = del(box)
	return box
end


function Dg(func,x,period)#ヤコビアン
	dimension = 2
    return strictJ2(func,x,period) - Matrix{BigFloat}(I,dimension,dimension)
end

function g(func,x,period)#
	fx = func(x)
	for i = 1:period-1
		fx = func(fx)
	end
    return  fx - x
end

function del(points,point)#まためんどくさいものを…取ってきた時にやれば良いのでは．　
	f(x) = mapslices(norm,x,dims = 1)　　　
		println(points .- point)　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
		boolean = f(points .- point).> 10^(-10)
		println(typeof(boolean))
		boolean = [boolean[i] for i =  1 : length(boolean)]
		println(typeof(boolean))
		points = points[:,boolean]
		println(points)
		points = hcat(points,point)
		println(points)
	return points  
end

#function del2(points,point)#こっちの方が早いかも.

#end

function G_Hat(matrix)#周期点に対応するヤコビアンに対応する,全てのSの候補を返す、つもり.思いつかないので2 * 2の場合ってことで許して.一応ベキ乗法のも後で組んでみることにする.
	mat2 =  convert(Array{Float64,dimension},matrix)
	w,evec = LinearAlgebra.eigen(mat2)
	Lambda = LinearAlgebra.Diagonal(w)#diagon
	eye = Matrix{BigFloat}(I,dimension,dimension)#単位行列
	S = eye
	reserve = evec * (S * Diagonal(w) - eye) * inv(evec)
	println(reserve)
	for i = 1 : dimension
		println(dimension)
		println(reserve)
		if abs(w[i]) > 1
			S[i,i] = -S[i,i]
			G = evec * (S * Diagonal(w) - eye) * inv(evec)
			println("s")
			reserve = vcat(reserve,G)#あらかじめS * Lambda - In を返すようにした.
			println(reserve)
			for j = i+1 : dimension
				if abs(w[j]) > 1
					S[j,j]  = -S[j,j]
					G = evec * (S * Diagonal(w) - eye) * inv(evec)
					reserve = vcat(reserve,G)#あらかじめS * Lambda - In を返すようにした.
					S[i,i] = -S[i,i]
					G = evec * (S * Diagonal(w) - eye) * inv(evec)
					reserve = vcat(reserve,G)#あらかじめS * Lambda - In を返すようにした.
				end
			end
		end
	#	G = evec * (S * Diagonal(w) - eye) * inv(evec)
	#	reserve = vcat(reserve,G)#あらかじめS * Lambda - In を返すようにした.
	end
	println(reserve)
	return reserve
end


function polardecomposition(G_hat)#STを求める際に極分解を行う.
	 F = svd(G_hat)
	 Q = F.U * F.Vt
	 println(Q * Q')
	 B = F.V * Diagonal(F.S) * F.Vt
	 return Q #帰ってくるのはQ.後でCに直す.
end

function Cs(Q)#極分解を行ったものの中からpolardecompo
	 C = - Q'
	return C
end

function DLMethod(func,x,period,C,beta)#Cを求めたら次はDL法を用いて周期点を見つける.xはp-1周期点とか(p+1)とか初期の点
	for  i =  1 : 10000
		print(x)
		x = x  + inv(beta * norm(g(func,x,period)) * C' - Dg(func,x,period)) * g(func,x,period)
		print(norm(g(func,x,period)) )
		println((1/beta) )
		if norm(g(func,x,period)) < (1/(beta))* 10#一桁上まで
			return x
		end
		if norm(g(func,x,period)) > 1
			break
		end
	end
	println(norm(g(func,x,period)) )
	return 2019
end

#function refine(beta,cmap.U,)#DL法で求まった方法をさらに精度をあげる
#end

function pointgetter(box,period,point)#周期点一つからたくさん求める
	for i  = 1:period-1
		point = cmap.U(point)
		box = hcat(box,point)
	end
	return box
end 
##########################################
function otherperiod(period,cmap,seeds)#一つ上の周期軌道を求めようとするプログラム.mainをサブルーチン化した.
	for i = 2 : Int(length(seeds[1,:]))#点ごとにstablized matrix.を求めてx_iの収束先を求めてp+1周期の点を手に入れる.あまりに計算時間が長かったら捨ててもいい.
		 x = seeds[1:2,i]
		 println(x)
		 matrix =  strictJ2(cmap.U,x,period)#p-1周期点におけるDg(x)
		 st = G_Hat(st)#stability matrixの元を返す.
		 print(st)
		 for j = 1 : Int(length(st)/(dimension * dimension))#reserveに入っているG_Hat一つ一つに対応してCsを得る.その後DLMethodで発散しないものから4周期点
		 	a = st[2 * j-1 : 2 * j,:]
		 	Q = polardecomposition(a)
		 	C  = Cs(Q)#一つのCsを得る.
		 	beta = 10^4#ここの値を低い周期で実験して調節する
		 	#if norm(g(cmap.U,x,period)) < 1/beta
		 	#	beta  = beta
		 	#end
		 	point = DLMethod(cmap.U,x,period+1,C,beta)#さらにこのCsについてsequenceを得る
		 	println(point)
		 	while true#どんどん精度を上げるところ
		 		if beta > 10^(16) || point == 2019
					println(point)
					break
				end	
		 		beta = beta * 2
		 		point = DLMethod(cmap.U,point,period+1,C,beta)#さらにこのCsについてsequenceを得る
		 		println(beta)
		 		println(point)
			end
			println(point)
			println(beta)
			if point != 2019
				points = hcat(points,point)
				println(points)
				println(typeof(points))
			end
		end
	end
	return points
end
#########################################

function main()
	setprecision = 20
	cmap = KScatteringMap
	boxs = randompp(4,cmap.U)#初期設定の部分.のちに追加する場合は、 
	println(boxs[:,1])
	boxs = del(boxs,boxs[:,2])
	period = 4
	points = [-1.2,0.0]
	otherperiod(period,cmap,boxs)#サブルーチン化したものを試す場所

	for i = 2 : Int(length(boxs[1,:]))#点ごとにstablized matrix.を求めて、sequenceを出して、p+1周期の点を手に入れる.
		 x = boxs[1:2,i]
		 println(x)
		 matrix =  strictJ2(cmap.U,x,period)#p-1周期点におけるDg(x)
		 st = G_Hat(matrix)#stability matrixの元を返す.
		 print(st)
		 #otherperiod()
		 for j = 1 : Int(length(st)/(dimension * dimension))#reserveに入っているG_Hat一つ一つに対応してCsを得る.その後DLMethodで発散しないものから4周期点
		 	a = st[2 * j-1 : 2 * j,:]
		 	Q = polardecomposition(a)
		 	C  = Cs(Q)#一つのCsを得る.
		 	beta = 10^4
		 	#if norm(g(cmap.U,x,period)) < 1/beta
		 	#	beta  = beta
		 	#end
		 	point = DLMethod(cmap.U,x,period+1,C,beta)#さらにこのCsについてsequenceを得る
		 	println(point)
		 	while true#どんどん精度を上げるところ
		 		if beta > 10^(16) || point == 2019
					println(point)
					break
				end	
		 		beta = beta * 2
		 		point = DLMethod(cmap.U,point,period+1,C,beta)#さらにこのCsについてsequenceを得る
		 		println(beta)
		 		println(point)
			end
			println(point)
			println(beta)
			if point != 2019
				points = hcat(points,point)
				println(points)
				println(typeof(points))
			end
		end
	end
	#ここまでで初期設定が終了する。
	print(points)
	boxs = hcat(boxs,points)
	points = box[1:2,:]
		#次に求めた点を用いて、次の周期点を求める。

end

main()
