//贝塞尔曲线绘制(javascript)
//Version: 1.0
//License: GPLv3





{
	/*
	//逻辑演算如下：=
	Zr bezierFunc {
		Rr (father) rec (Q)
		Rr (child) rec (Q)
		Rr (Q) iter (Dlm,L) dev (Mem)
		Rr (AfterL) rec (Li) dev (Lif)
	}
	iter (T) dev (R) rec (Y) dev (M)
	 */

	//计算阶乘
	let 阶乘 = ( n ) => {
		let res = 1;
		for( let i = 1; i <= n; i++ ){
			res *= i;
		}
		return res;
	}
	//一个总的绘制函数
	let beisaier = (obj) => {
		//定义一个绘制函数,point=[[xyz1],[xyz2],[xyz3]...]
		let draw = (point) => {
			let ps = []
			//先计算这些控制点构成的曲线长度,然后采样的点数=Math.round(长度-0.5)
			//计算贝塞尔曲线的1阶导数,并且是在0和1两个位置
			let dfdx = (t) => {
				//根据公式返回
				//C^(k)(t) = n(n-1)...(n-k+1)∑[i=0,n-k][B[i,n-k]p[i]^(k)],p[i]^(k)=p[i+1]^(k-1)-p[i]^(k-1)
				//B[i,n]=(n!)/(i!(n-i)!)*(1-t)^(n-i)*t^i
				//现在只需要得出C^(1)(t)的解析式,
				let dx = 0
				let dy = 0
				let dz = 0
				let n = point.length
				for(let i=0;i<n-1;i++){
					let p = point[i]
					let prex = n*(point[i+1][0]-p[0])
					let prey = n*(point[i+1][1]-p[1])
					let prez = n*(point[i+1][2]-p[2])
					let pre = 阶乘(n-1)/(阶乘(i)*阶乘(n-i-1))*Math.pow(1-t,n-i-1)*Math.pow(t,i)
					dx += pre*prex
					dy += pre*prey
					dz += pre*prez
				}
				return [dx,dy,dz]
			}
			//三阶导数
			let d2fdx3 = (t) => {
				//.....
				let dx = 0
				let dy = 0
				let dz = 0
				let n = point.length
				for(let i=0;i<n-2;i++) {
					let p = point[i]
					let prex = (n-1)*n*(point[i+2][0]-2*point[i+1][0]+p[0])
					let prey = (n-1)*n*(point[i+2][1]-2*point[i+1][1]+p[1])
					let prez = (n-1)*n*(point[i+2][2]-2*point[i+1][2]+p[2])
					let pre = 阶乘(n-1)/(阶乘(i)*阶乘(n-i-1))*Math.pow(1-t,n-i-1)*Math.pow(t,i)
					dx += pre*prex
					dy += pre*prey
					dz += pre*prez
				}
				return [dx,dy,dz]
			}
			//魔改这个函数,使得曲线的d2fdx3的解析式可以更好的求出来,这个函数的解析式是: d2fdx3(t) = dfdx(t) + dfdx(t)^2*d2fdx3(t) + dfdx(t)^3*d3fdx3(t) + ... = dfdx(t) + dfdx(t)^2*d2fdx3(t) + dfdx(t)^3*d3fdx3(t) + ...


			//直接定积分求长度,
			let L = (t) => {
				//L = ∫[0,t]|dfdx(t)|dt = C(t) - C(0)
				//又 C(t) = ∑[i=0,n]B[i,n](t)P[i]
				let dx = 0
				let dy = 0
				let dz = 0
				let n = point.length
				let a = 0
				let b = t
				let mida = (2*a+b)/3
				let midb = (a+2*b)/3
				let d1 = dfdx(a)
				d1 = Math.sqrt(Math.pow(d1[0],2)+Math.pow(d1[1],2)+Math.pow(d1[2],2))
				let d2 = dfdx(mida)
				d2 = Math.sqrt(Math.pow(d2[0],2)+Math.pow(d2[1],2)+Math.pow(d2[2],2))
				let d3 = dfdx(midb)
				d3 = Math.sqrt(Math.pow(d3[0],2)+Math.pow(d3[1],2)+Math.pow(d3[2],2))
				let d4 = dfdx(b)
				d4 = Math.sqrt(Math.pow(d4[0],2)+Math.pow(d4[1],2)+Math.pow(d4[2],2))
				//使用魔法进行数值积分计算
				//上面的曲线可以通过一个特殊的公式求解,因为这个曲线没有一个点的曲率,所以可以求出曲线的长度
				return (d1+d2*3+d3*3+d4)*(b-a)/8
			}
			//计算t对应的点
			let t2point = (t) => {
				let x = 0
				let y = 0
				let z = 0
				let n = point.length
				for(let i=0;i<n;i++){
					let p = point[i]
					let px = p[0]
					let py = p[1]
					let pz = p[2]
					let pre = 阶乘(n)/(阶乘(i)*阶乘(n-i))*Math.pow(1-t,n-i)*Math.pow(t,i)
					x += pre*px
					y += pre*py
					z += pre*pz
				}
				return [x,y,z]
			}
			//计算nextT
			let nextT = (t,s) => {
				//会算出dx/dy/dz最小的t值
				let c1 = dfdx(t)
				let d = L(t) - s
				let c = t2point(t)
				//console.log([d,s])
				let g1 = [
					2*d*Math.abs(c1[0]),
					2*d*Math.abs(c1[1]),
					2*d*Math.abs(c1[2])
				]
				let mx = 1/g1[0]
				let my = 1/g1[1]
				let mz = 1/g1[2]
				//console.log(t)
				return Math.min(mx,my,mz)+t
			}
			let maxL = L(1)
			let lp
			//console.log(maxL)
			//令s=0,牛顿迭代直到s>=L,
			let t = 0.1
			for(let i = 0;i<maxL;i++){
				let p = t2point(t)
				p = [
					Math.round(p[0]),
					Math.round(p[1]),
					Math.round(p[2])
				]
				if(t>1) t=1
				if (lp) {
					//如果是相等的,则不push
					if (p[0]!=lp[0]||p[1]!=lp[1]||p[2]!=lp[2]) {
						lp = p
						ps.push(p)
					}
				} else {
					lp = p
					ps.push(p)
				}
				t = nextT(t,i+1)
				//console.log(t)
			}
			return ps
		}
		let father = draw(obj.father)
		let child = draw(obj.child)
		let aim = obj.aim || null
		if(aim==null){
			let minx = 0
			let miny = 0
			let minz = 0
			let maxx = 0
			let maxy = 0
			let maxz = 0
			for(let i=0;i<child.length;i++) {
				let p = child[i]
				if(p[0]<minx) minx = p[0]
				if(p[1]<miny) miny = p[1]
				if(p[2]<minz) minz = p[2]
				if(p[0]>maxx) maxx = p[0]
				if(p[1]>maxy) maxy = p[1]
				if(p[2]>maxz) maxz = p[2]
			}
			aim = [
				(minx+maxx)/2,
				(miny+maxy)/2,
				(minz+maxz)/2
			]
		}
		function b () {
			return obj.block[Math.round(Math.random()*(obj.block.length-0.5))]
		}
		//以father的轨迹,绕着aim旋转child
		let n = father.length
		let ps = []
		let last
		for(let i=0;i<n;i++){
			for(let j=0;j<child.length;j++){
				let fx = father[i][0]
				let fy = father[i][1]
				let fz = father[i][2]
				let ex = aim[0]
				let ey = aim[1]
				let ez = aim[2]
				let cx = child[j][0]
				let cy = child[j][1]
				let cz = child[j][2]
				let x = Math.round(ex-cx+fx)
				let y = Math.round(ey-cy+fy)
				let z = Math.round(ez-cz+fz)
				//如果xyz又NaN,则跳过
				if(isNaN(x)||isNaN(y)||isNaN(z)) continue
				let bs = b()
				if (!last) {
					last = [x,y,z,bs]
					ps.push(last)
				} else {
					if (last[0]!=x||last[1]!=y||last[2]!=z) {
						last = [x,y,z,bs]
						ps.push(last)
					}
				}
			}
		}
		return ps
	}
	/*let res = beisaier({
		father:[[-41,51,67],[-16,27,96],[4,4,4],[70,76,45],[78,26,117]],
		child:[[73,36,-84],[34,45,14],[26,25,22],[35,46,72]],
		aim:null,
		block: ['grass 0','wool 12','log 2','stone 0']
	})*/

	exports.bezier = beisaier
	setTimeout(()=>{
		console.log('bezier.js loaded')
	})
}

