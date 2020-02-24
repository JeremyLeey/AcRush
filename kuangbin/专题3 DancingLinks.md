[toc]

# 精准覆盖问题
中文直译跳舞链,用于解决精确覆盖问题,数学模型为
```
在一个全集X中若干子集的集合为S,精确覆盖是指,S的子集S*,满足X中的每一个元素在S*中恰好出现一次,即
1) S*中任意两个集合没有交集,即X中的元素在S*中出现最多一次
2) S*中集合的全集为X,即X中的元素在S*中出现最少一次
合二为一,即X中的元素在S*中出现恰好一次
```
解决的思路是dfs,Dancing Links的含义在于使用交叉十字循环双向链这样一种数据结构来比较方便的实现列和行的插入与删除

参考文献:
1. https://www.cnblogs.com/grenet/p/3145800.html

# 20.2.12 hihoCoder 1317 精准覆盖模板
```
struct DLX {
	const int head = 0;
	int n, m;
	// U[i] 第i个节点的up指针
	int U[mm], D[mm], L[mm], R[mm];
	// Row[i],Col[i] 第i个节点所对应的行号,列号
	int Row[mm], Col[mm];
	// 辅助数组,H[i]代表第i行的第一个元素节点,colcnt[i]第i列元素个数
	int H[maxn], colcnt[maxn];
	int size; // 节点个数
	void init(int n_, int m_) {
		// 初始化表头(第一行)
		n = n_; m = m_;
		for (int i = 0; i <= m; i++) {
		    // 上下指针暂时没用到
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
		}
		R[m] = 0, L[0] = m;
		for (int i = 1; i <= n; i++) {
			H[i] = -1;
		}

		size = m + 1;
	}

	// 插入点(r, c)
	void link(int r, int c) {
		Col[size] = c;
		Row[size] = r;
		
		U[size] = U[c]; D[U[c]] = size;
		D[size] = c; U[c] = size;

		if (H[r] < 0) {
			// 这一行没有元素
			H[r] = L[size] = R[size] = size;
		}
		else {
			L[size] = L[H[r]]; R[L[H[r]]] = size;
			L[H[r]] = size; R[size] = H[r];
		}
		size++;
	}
    // 删除第col列
	void remove(int col) {
		R[L[col]] = R[col];
		L[R[col]] = L[col];
		// 删除该列元素所在的行的元素
		for (int i = D[col]; i != col; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = U[j]; D[U[j]] = D[j];
				colcnt[C[j]]--;
			}
		}
	}
    // 恢复第col列
	void resume(int col) {
		R[L[col]] = col, L[R[col]] = col;
		for (int i = D[col]; i != col; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = j; D[U[j]] = j;
				coln[C[j]]++;
			}
		}
	}
    // 搜索过程, k代表搜索层数
	int dance(int k) {
	    // 没有要覆盖的列
		if (R[head] == head) {
			return 1;
		}
		
		/** 这个地方有一个优化,每次选择包含1最小的列
		    因为这样接下来的循环就会次数少
		    这也是为什么上面存储了colcnt数组
		    int minn=maxn,minpos=-1;
            for(int i=R[0];i!=0;i=R[i]){
                if(minn>coln[i]){
                    minn=coln[i];  minpos=i;
                }
        **/
        
		int C = R[head];
		remove(C);
		for (int i = D[C]; i != C; i = D[i]) {
			// 枚举包含该列的每一行
			for (int j = R[i]; j != i; j = R[j]) {
				remove(Col[j]);
			}
			if (dance(k + 1)) return 1;
			// 回溯恢复,注意这里要逆向
			for (int j = L[i]; j != i; j = L[j]) {
				resume(Col[j]);
			}
		}
		resume(C);
		return 0;
	}
}dlx;
```

# 20.2.13 ZOJ 3209 TreasureMap 精准覆盖模板题
给定一些矩形方块,同时指定一块矩形区域,问你能否从这些方块中选出一些,使它们拼在一起恰好覆盖这个矩形区域,且没有重叠。如果可以,返回最少的方块数,否则返回-1

---

这是一道精准覆盖的问题,对其进行数学抽象,可以把每一个方块看成一行,把指定区域划分为1x1的小方块,然后把每一个小方块看成一列,这样就成了dancing link的模板题;一开始TLE了一发,后来加了每次选择包含元素最少的那一列这一优化就可以通过了

```
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

const int M = 550;
const int N = 1000;
const int head = 0;
int t, ans;
struct DLX {
	int U[M*N], D[M*N], L[M*N], R[M*N];
	int Row[M*N], Col[M*N];
	int H[M], S[N];
	int sz, m, n;
	void init(int m_, int n_) {
		// 初始化表头
		m = m_; n = n_;
		for (int i = 0; i <= n; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = n; R[n] = 0;
		sz = n + 1;
		for (int i = 1; i <= m; i++) {
			H[i] = -1;
		}
	}

	void add(int r, int c) {
		// r行c列添加一个元素
		Row[sz] = r, Col[sz] = c;
		S[c]++;

		U[sz] = U[c]; D[U[c]] = sz;
		D[sz] = c; U[c] = sz;
		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			L[sz] = L[H[r]]; R[L[H[r]]] = sz;
			L[H[r]] = sz; R[sz] = H[r];
		}
		sz++;
	}

	void remove(int c) {
		// 删除第col列元素
		R[L[c]] = R[c]; L[R[c]] = L[c];
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				D[U[j]] = D[j];
				U[D[j]] = U[j];
				S[Col[j]]--;
			}
		}
	}

	void resume(int c) {
		R[L[c]] = c; L[R[c]] = c;
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				D[U[j]] = j;
				U[D[j]] = j;
				S[Col[j]]++;
			}
		}
	}

	void dance(int depth) {
		if (depth > ans) {
			return;
		}
		if (R[head] == head) {
			ans = min(ans, depth);
		}
		// 选择优化(一定要加)
		int C = R[head];
		for (int i = R[head]; i != head; i = R[i]) {
			if (S[i] < S[C]) {
				C = i;
			}
		}

		remove(C);
		for (int i = D[C]; i != C; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				remove(Col[j]);
			}
			dance(depth+1);
			for (int j = L[i]; j != i; j = L[j]) {
				resume(Col[j]);
			}
		}
		resume(C);
	}
	
}dlx;


int main() {
	// freopen("E://input.txt", "r", stdin);
	scanf("%d", &t);
	while (t--) {
		int n, m, p;
		scanf("%d %d %d", &n, &m, &p);
		dlx.init(p, m*n);
		for (int i = 1; i <= p; i++) {
			int x1, y1, x2, y2;
			scanf("%d %d %d %d", &x1, &y1, &x2, &y2);
			for (int y = y1; y < y2; y++) {
				for (int x = x1; x < x2; x++) {
					dlx.add(i, n*y+x+1);
				}
			}
		}
		ans = 0x3f3f3f3f;
		dlx.dance(0);
		printf("%d\n", ans < 0x3f3f3f3f ? ans : -1);
	}
	return 0;
}
```

# 重复覆盖问题
所谓的重复覆盖是指
```
在一个全集X中若干子集的集合为S,重复覆盖是指,S的子集S*,满足X中的每一个元素在S*中至少出现1次
```
由此可以看出,精准覆盖是重复覆盖的特殊情况。重复覆盖同样使用dancing link进行建模，模板中init函数与addlink函数一致,区别在于remove,resume以及dance函数。
```
思考精准覆盖的删除过程,当选中一列时,由于该列要求有且仅有1个1,因此需要删除该列包括的所有行,同时当枚举每一行时,要删除这一行包括的列。

由于重复覆盖该列要求至少1个，因此不需要删除该列对应的所有行，只能删除该列
```
由于搜索过程中矩阵规模的下降程度会减慢，因此需要加入强剪枝操作，即采用迭代加深搜索，利用估价函数f提前进行剪枝判断。

# 20.2.14 HDU 3498 重复覆盖模板题
地图上有n个小兵的位置,并且告诉你一些小兵之间的相邻关系，玩家攻击一名小兵时会同样攻击这名小兵的相邻单位，问你最少需要攻击多少次，能将地图上的小兵全部杀死。

---

列出小兵的相邻矩阵,m[i][j]为1表示第i个小兵和第j个小兵相邻，那么问题转化为了从中选择最少的行，使其的每一列至少有1个1（即至少被攻击一次），这是一道重复覆盖的模板题

```
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

const int oo = 0x3f3f3f3f;
int N, M, ans, board[60][60];

struct DLX {
	static const int maxN = 60;
	static const int maxM = 60;
	static const int maxMN = maxM * maxN;
	static const int head = 0;
	int n, m; // n行m列
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxN], S[maxM];
	int sz;

	void init(int n_, int m_) {
		m = m_;
		n = n_;
		for (int i = 0; i <= m; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
		}
		L[0] = m; R[m] = 0;
		for (int i = 1; i <= n; i++) {
			H[i] = -1;
		}
		sz = m + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r; Col[sz] = c;
		S[c]++;

		D[U[c]] = sz; U[sz] = U[c];
		D[sz] = c; U[c] = sz;
		if (H[r] == -1) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz; L[sz] = L[H[r]];
			R[sz] = H[r]; L[H[r]] = sz;
		}
		sz++;
	}

    // 重复覆盖,c为元素编号而非列号
    // 即删除该元素所在的列 包括它自身和列头指针
	void repeat_remove(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = L[i], R[L[i]] = R[i];
		}
	}
	// ps:其实这里有一个困惑的地方：并没有删除元素编号为c的这个节点的左右关系，但是并不影响

	void repeat_resume(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = R[L[i]] = i;
		}
	}

    // A*的估价函数
	int f() {
		/*
		 这个剪枝的思想是A*搜索中的估价函数。当前搜索的深度代表已经选择的行数，即对于当前递归深度
		 K下的矩阵，估计其最少还需要多少部才能出解，如果已经超过了设置的最大搜索深度，则直接返回
		 这个估算思路是 如果将能够覆盖当前列的所有行全部选中，去掉这些行能够覆盖到的列，将这个操作
		 作为一层深度，重复此操作直到所有列全部出解的深度是多少。
		*/
		bool vv[maxM];
		int ret = 0;
		for (int c = R[head]; c != head; c = R[c]) {
			vv[c] = true;
		}
		for (int c = R[head]; c != head; c = R[c]) {
			if (vv[c]) {
				ret++;
				vv[c] = false;
				for (int i = D[c]; i != c; i = D[i]) {
					for (int j = R[i]; j != i; j = R[j]) {
						vv[Col[j]] = 0;
					}
				}
			}
		}
		return ret;
	}

	bool repeat_dance(int depth, int maxdepth) {
		if (depth + f() > maxdepth) {
			return false;
		}
		if (R[head] == head) {
			ans = min(ans, depth);
			return true;
		}
		int C = R[head];
		for (int i = R[head]; i != head; i = R[i]) {
			if (S[i] < S[C]) {
				C = i;
			}
		}
        
        // 枚举该列下的每一行
		for (int i = D[C]; i != C; i = D[i]) {
		    // 先删除该列
			repeat_remove(i);
			// 删除该行包括的其余列
			for (int j = R[i]; j != i; j = R[j]) {
				repeat_remove(j);
			}
			if (repeat_dance(depth + 1, maxdepth)) return true;
			for (int j = L[i]; j != i; j = L[j]) {
				repeat_resume(j);
			}
			repeat_resume(i);
		}
		return false;

	}
}dlx;

int main() {
	// freopen("E://input.txt", "r", stdin);
	while (~scanf("%d %d", &N, &M)) {
		memset(board, 0, sizeof board);
		ans = oo;
		dlx.init(N, N);
		for (int i = 0; i < M; i++) {
			int u, v;
			scanf("%d %d", &u, &v);
			board[u][v] = board[v][u] = 1;
		}
		for (int i = 1; i <= N; i++) {
			for (int j = 1; j <= N; j++) {
				if (i == j || board[i][j] == 1) {
					dlx.addlink(i, j);
				}
			}
		}

		for (int maxd = 0; maxd <= N; maxd++) {
			if (dlx.repeat_dance(0, maxd)) {
				printf("%d\n", maxd);
				break;
			}
		}
	}
	return 0;
}

```

# 20.2.14 HDU 2295 Radar 2分+重复覆盖IDA* (★★★)
给你n个城市和m个雷达的坐标，这m个雷达都拥有同样的半径，求一个最小的半径，使每个城市都能被不超过k个雷达覆盖。

---

一种很常见的思路，对要求解的变量范围进行二分。当半径为r时，能够得到一个雷达和城市的覆盖关系图（01矩阵），然后问题就转化了对这个矩阵进行一个重复覆盖问题的求解。

```
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

struct DLX {
	static const int maxN = 60;
	static const int maxM = 60;
	static const int maxMN = maxM * maxN;
	static const int head = 0;
	int n, m; // n行m列
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxN], S[maxM];
	int sz;

	void init(int n_, int m_) {
		m = m_;
		n = n_;
		for (int i = 0; i <= m; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
		}
		L[0] = m; R[m] = 0;
		for (int i = 1; i <= n; i++) {
			H[i] = -1;
		}
		sz = m + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r; Col[sz] = c;
		S[c]++;

		D[U[c]] = sz; U[sz] = U[c];
		D[sz] = c; U[c] = sz;
		if (H[r] == -1) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz; L[sz] = L[H[r]];
			R[sz] = H[r]; L[H[r]] = sz;
		}
		sz++;
	}

	void repeat_remove(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = L[i], R[L[i]] = R[i];
		}
		//L[R[c]] = L[c]; R[L[c]] = R[c];
	}

	void repeat_resume(int c) {
		//L[R[c]] = R[L[c]] = c;
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = R[L[i]] = i;
		}
	}

	int f() {
		bool vv[maxM];
		int ret = 0;
		for (int c = R[head]; c != head; c = R[c]) {
			vv[c] = true;
		}
		for (int c = R[head]; c != head; c = R[c]) {
			if (vv[c]) {
				ret++;
				vv[c] = false;
				for (int i = D[c]; i != c; i = D[i]) {
					for (int j = R[i]; j != i; j = R[j]) {
						vv[Col[j]] = 0;
					}
				}
			}
		}
		return ret;
	}

	bool repeat_dance(int depth, int maxdepth) {
		if (depth + f() > maxdepth) {
			return false;
		}
		if (R[head] == head) {
			return true;
		}
		int C = R[head], temp = S[C];
		for (int i = R[head]; i != head; i = R[i]) {
			if (S[i] < temp) {
				temp = S[i];
				C = i;
			}
		}

		//
		for (int i = D[C]; i != C; i = D[i]) {
			repeat_remove(i);
			for (int j = R[i]; j != i; j = R[j]) {
				repeat_remove(j);
			}
			if (repeat_dance(depth + 1, maxdepth)) return true;
			for (int j = L[i]; j != i; j = L[j]) {
				repeat_resume(j);
			}
			repeat_resume(i);
		}
		return false;

	}
}dlx;

const double delta = 1e-7;

int cityX[55], cityY[55];
int radarX[55], radarY[55];

double distance(double x1, double y1, double x2, double y2) {
	return (x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2);
}

int main() {
	// freopen("E://input.txt", "r", stdin);
	int t;
	scanf("%d", &t);
	while (t--) {
		int n, m, k;
		scanf("%d %d %d", &n, &m, &k);
		for (int i = 1; i <= n; i++) {
			scanf("%d %d", &cityX[i], &cityY[i]);
		}
		for (int i = 1; i <= m; i++) {
			scanf("%d %d", &radarX[i], &radarY[i]);
		}
		double l = 0, h = 2000, mid = 0;
		while (h - l > delta) {
			mid = (l + h) / 2;
			// 雷达是行 城市是列
			dlx.init(m, n);
			for (int i = 1; i <= m; i++) {
				for (int j = 1; j <= n; j++) {
					if (distance(radarX[i], radarY[i], cityX[j], cityY[j]) <= mid*mid) {
						dlx.addlink(i, j);
					}
				}
			}
			if (dlx.repeat_dance(0, k)) {
				h = mid;
			}
			else {
				l = mid;
			}
		}
		printf("%.6llf\n", mid);
	}
	return 0;
}
```

# 20.2.15 POJ 1084 Square Destroyer 建模+重复覆盖IDA*
给定一个nxn的网格,网格的每个单位都是由火柴拼成的,并且从上到下从左到右火柴依次标号。现在从中拿掉一些火柴，将会破怪一些正方形，现在问你还需要最少再拿几根，可以破坏全部的正方形。

---

因为是Dancing Link专题，自然就会往Dancing Link去想，模型也很好抽象,即火柴作为行，每一个正方形作为列,当去掉这根火柴该正方形被破坏,则值为1否则值为0。

但是这题麻烦的地方有两个，一个是矩阵的建立，另一个是题目会提前告诉你拿掉一些火柴。矩阵的建立采用循环每一个正方形(先循环边长再循环左上角坐标),然后由左上角坐标确定组成边的火柴编号,然后依次循环四条边,将找到的每一条火柴的对应位置1即可。第二个就是题目会告你提前拿掉几根火柴，一旦拿掉一根火柴，那么该火柴影响的矩阵则全部不需要判断，即该行包含的列全部无效。这里采用一个标记(代表该列已经被删除)，然后再添加节点的时候判断该标记是否有效即可。注意需要修改表头的左右指针以表示该列被完全删除。

```
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <vector>

using namespace std;


const int total[][2] = {0,0, 4,1, 12,5, 24,14, 40,30, 60,55}; // n=1,2,3,4,5的火柴总数(0)与正方形总数(1)
const int maxM = 65, maxN = 60;
const int maxMN = maxM* maxN;

int n, m, k,board[maxM][maxN]; // board[i][j], 第i根火柴影响第j个正方形
bool isDeleted[maxN];

struct DLX {
	int M, N, sz;
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxM], S[maxN];
	void init(int row, int col) {
		M = row, N = col;
		for (int i = 0; i <= N; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = N, R[N] = 0;
		memset(H, -1, sizeof H);
		sz = N + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r, Col[sz] = c;
		S[c]++;
		D[U[c]] = sz, U[sz] = U[c];
		U[c] = sz, D[sz] = c;

		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz, L[sz] = L[H[r]];
			R[sz] = H[r], L[H[r]] = sz;
		}

		sz++;
	}

	void repeat_remove(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			R[L[i]] = R[i];
			L[R[i]] = L[i];
		}
	}

	void repeat_resume(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			R[L[i]] = L[R[i]] = i;
		}
	}

	int f() {
		bool tag[maxN];
		int ret = 0;
		for (int i = 0; i <= N; i++) tag[i] = true;
		for (int c = R[0]; c != 0; c = R[c]) {
			if (tag[c]) {
				ret++;
				tag[c] = false;
				for (int i = D[c]; i != c; i = D[i]) {
					for (int j = R[i]; j != i; j = R[j]) {
						tag[Col[j]] = false;
					}
				}
			}
		}
		return ret;
	}

	bool repeat_dance(int depth, int maxdepth) {
		if (depth + f() > maxdepth) return false;
		if (R[0] == 0) {
			return true;
		}
		int c = R[0];
		for (int i = R[0]; i != 0; i = R[i]) {
			if (S[i] < S[c]) {
				c = i;
			}
		}

		for (int i = D[c]; i != c; i = D[i]) {
			repeat_remove(i);
			for (int j = R[i]; j != i; j = R[j]) {
				repeat_remove(j);
			}
			if (repeat_dance(depth + 1, maxdepth)) return true;
			for (int j = L[i]; j != i; j = L[j]) {
				repeat_resume(j);
			}
			repeat_resume(i);
		}
		return false;
	}
}dlx;

void sol(int x, int y, int l, int sz) {
	// 以(i, j)为左上顶点,l为边长的矩形 计算每条火柴的编号即可
	int x0 = (2 * n + 1)*x + 1 + y, t = x0; // x0: x行的第一根火柴编号
	
	// 上边的火柴编号
	for (int i = 1; i <= l; i++) {
		board[t][sz] = 1;
		t++;
	}
	// 下边的火柴编号
	t = (2 * n + 1)*l + x0;
	for (int i = 1; i <= l; i++) {
		board[t][sz] = 1;
		t++;
	}
	// 左边的火柴编号
	t = x0 + n;
	for (int i = 1; i <= l; i++) {
		board[t][sz] = 1;
		t += 2*n + 1;
	}

	// 右边的火柴编号
	t = x0 + n + l;
	for (int i = 1; i <= l; i++) {
		board[t][sz] = 1;
		t += 2*n + 1;
	}

}

void preProcess() {
	int sz = 1; // 矩形编号
	for (int l = 1; l <= n; l++) {
		for (int i = 0; i + l <= n; i++) {
			for (int j = 0; j + l <= n; j++) {
				// 以(i, j)为左上顶点,l为边长的矩形
				sol(i, j, l, sz++);
			}
		}
	}
}

int main() {
	// freopen("E://input.txt", "r", stdin);
	int t;
	scanf("%d", &t);
	while (t--) {
		memset(board, 0, sizeof board);
		scanf("%d %d", &n, &k);
		
		preProcess();
		dlx.init(total[n][0], total[n][1]);

		memset(isDeleted, 0, sizeof isDeleted);
		for (int i = 0; i < k; i++) {
			int x;
			scanf("%d", &x);
			for (int j = 1; j <= total[n][1]; j++) {
				if (board[x][j] == 1 && isDeleted[j] == 0) {
				    // 标记该列已经被删除且修改表头的左右指针
					isDeleted[j] = 1;
					dlx.L[dlx.R[j]] = dlx.L[j]; dlx.R[dlx.L[j]] = dlx.R[j];
				}
			}
		}
		// 在插入节点的时候判断是否有效
		for (int i = 1; i <= total[n][0]; i++) {
			for (int j = 1; j <= total[n][1]; j++) {
				if (!isDeleted[j] && board[i][j] == 1) {
					dlx.addlink(i, j);
				}
			}
		}
		for (int maxd = 0; ;maxd++) {
			if (dlx.repeat_dance(0, maxd)) {
				printf("%d\n", maxd);
				break;
			}
		}
	}
	return 0;
}
```

# 20.2.15 POJ 3740 Easy Finding 精确覆盖纯裸题
# 20.2.16 POJ 3074 Sudoku DancingLink解数独
行代表问题的所有情况,列代表约束条件。以九宫格数独为例，一共有81个格子,每个格子取值为1~9,所以一共有729种情况。即DancingLinks的有729行。列分为四种

```
1)[1, 81]列分别对应了81个格子是否被放置了数字(1为放置,0未放置)
2)[82, 81*2]列分别对应了9行每行[1,9]个数字的放置情况
3)[163, 81*3]列分别对应了9列每列[1,9]个数字的放置情况
4)[244, 81*4]列分别对应了9个小宫格,每个小宫格[1,9]个数字的放置情况
```
DancingLink建立的过程即为顺序遍历数独,在对于(i, j)上有数字k只需插入一行,这行上有4列为1。对于没有填写的数字,则需要枚举1到9,把在(i, j)上填[1, 9]的情况都进行插入,这就转变成了一个精确覆盖问题。
```C++
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

const int maxM = 750, maxN = 330;
const int maxMN = maxM * maxN;
int board[maxM][maxN];
int d, ans[maxM];
int soduku[10][10];

struct DLX {
	static const int head = 0;
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxM], S[maxN];
	int m, n, sz;
	void init(int m_, int n_) {
		m = m_, n = n_;
		for (int i = 0; i <= n; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = n, R[n] = 0;
		memset(H, -1, sizeof H);
		sz = n + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r, Col[sz] = c;
		S[c]++;
		D[U[c]] = sz, U[sz] = U[c];
		D[sz] = c, U[c] = sz;

		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz, L[sz] = L[H[r]];
			R[sz] = H[r], L[H[r]] = sz;
		}
		sz++;
	}
	
	void remove(int c) {
		R[L[c]] = R[c], L[R[c]] = L[c];
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				D[U[j]] = D[j];
				U[D[j]] = U[j];
				S[Col[j]]--;
			}
		}
	}

	void resume(int c) {
		R[L[c]] = L[R[c]] = c;
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = D[U[j]] = j;
				S[Col[j]]++;
			}
		}
	}

	bool dance(int depth) {
		if (R[head] == head) {
			d = depth;
			return true;
		}
		int c = R[head];
		for (int i = R[head]; i != head; i = R[i]) {
			if (S[i] < S[c]) {
				c = i;
			}
		}

		remove(c);
		for (int i = D[c]; i != c; i = D[i]) {
			ans[depth] = Row[i];
			for (int j = R[i]; j != i; j = R[j]) {
				remove(Col[j]);
			}
			if (dance(depth + 1)) return true;
			for (int j = L[i]; j != i; j = L[j]) {
				resume(Col[j]);
			}
		}
		resume(c);
		return false;
	}

}dlx;
char seq[90];

// 在该行的四个列上分别置为1
void insert(int i, int j, int k) {
	// i 行号1~9, j列号1~9
	int row = 9 * (9 * (i - 1) + j - 1) + k;
	board[row][9 * (i - 1) + j] = board[row][81 + 9 * (i - 1) + k] = board[row][81 * 2 + 9 * (j - 1) + k] = 1;
	int cell = 3 * ((i - 1) / 3) + ((j - 1) / 3) + 1;
	board[row][81 * 3 + 9 * (cell-1) + k] = 1;
}

int main() {
	// freopen("E://input.txt", "r", stdin);
	while (~scanf("%s", seq)) {
		if (strcmp(seq, "end") == 0) break;
		memset(board, 0, sizeof board);
		memset(soduku, 0, sizeof soduku);
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				int idx = 9 * i + j;
				if (seq[idx] == '.') {
					for (int n = 1; n <= 9; n++) {
						insert(i+1, j+1, n);
					}
				}
				else {
					// (i, j)放置了一个数字
					soduku[i][j] = seq[idx] - '0';
					insert(i+1, j+1, seq[idx]-'0');
				}
			}
		}
		dlx.init(729, 324);
		for (int i = 1; i <= 729; i++) {
			for (int j = 1; j <= 81 * 4; j++) {
				if (board[i][j] == 1) {
					dlx.addlink(i, j);
				}
			}
		}

		if (dlx.dance(0)) {
			for (int i = 0; i < d; i++) {
				int rowNum = ans[i];
				int id = (rowNum - 1) / 9 + 1, val = rowNum - 9 * (id - 1);
				soduku[(id-1)/9][(id-1)%9] = val;
			}
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 9; j++) {
					printf("%d", soduku[i][j]);
				}
			}
			printf("\n");
		}
		else {
			printf("no answer\n");
		}
	}
	return 0;
}
```
用ans数组记录搜索得到的答案。则对于ans[i]代表的是选择第i行,第i行所代表的含义是就是第xx个格子取值为yy。把xx和yy计算出来即可。

# 20.2.16 FZU-1686 神龙的难题 重复覆盖
这题OJ崩了并不能提交,简单说一下DancingLink矩阵的建立。行是第i次攻击,列是怪物。第i次攻击可以杀死怪物j则board[i][j] = 1;枚举攻击区域的左上角,然后遍历攻击区域,对于区域内的怪物则标记为1。

# 20.2.17 ZOJ 3122 Sudoku 4阶数独
同POJ-3074,不过是4阶数独，处理的方式一样，不过这题有多组输入用例,且每个测试用例的输出有一个空行，注意最后一个测试用例后没有空行。
```
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;
const int maxM = 4100;
const int maxN = 1050;
const int maxMN = maxM * maxN;
char grid[20][20];
int board[maxM][maxN], ans[maxM];
int d;
void insert(int x, int y, int val) {
	int num = 16 * (x - 1) + y;
	int row = 16 * (num - 1) + val;
	int cell = 4 * ((x - 1) / 4) + ((y - 1) / 4) + 1;
	int col1 = num;
	int col2 = 256 + 16 * (x - 1) + val;
	int col3 = 512 + 16 * (y - 1) + val;
	int col4 = 768 + 16 * (cell - 1) + val;
	board[row][col1] = board[row][col2] = board[row][col3] = board[row][col4] = 1;
}

struct DLX {
	int m, n;
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxM], S[maxN];
	int sz;

	void init(int m_, int n_) {
		m = m_, n = n_;
		for (int i = 0; i <= n; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = n, R[n] = 0;
		memset(H, -1, sizeof H);
		sz = n + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r, Col[sz] = c;
		S[c]++;

		D[U[c]] = sz, U[sz] = U[c];
		U[c] = sz, D[sz] = c;

		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz, L[sz] = L[H[r]];
			R[sz] = H[r], L[H[r]] = sz;
		}
		sz++;
	}

	void remove(int c) {
		L[R[c]] = L[c], R[L[c]] = R[c];

		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = U[j], D[U[j]] = D[j];
				S[Col[j]]--;
			}
		}
	}

	void resume(int c) {
		L[R[c]] = R[L[c]] = c;
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = D[U[j]] = j;
				S[Col[j]]++;
			}
		}
	}

	bool dance(int depth) {
		if (R[0] == 0) {
			d = depth;
			return true;
		}
		int c = R[0];
		for (int i = R[0]; i != 0; i = R[i]) {
			if (S[i] < S[c]) {
				c = i;
			}
		}

		remove(c);
		for (int i = D[c]; i != c; i = D[i]) {
			ans[depth] = Row[i];
			for (int j = R[i]; j != i; j = R[j]) {
				remove(Col[j]);
			}
			if (dance(depth + 1)) return true;
			for (int j = L[i]; j != i; j = L[j]) {
				resume(Col[j]);
			}
		}
		resume(c);
		return false;
	}
}dlx;
char line[20];

void process(int i) {
	for (int j = 0; j < 16; j++) {
		if (line[j] == '-') {
			for (int val = 1; val <= 16; val++) {
				insert(i + 1, j + 1, val);
			}
		}
		else {
			insert(i + 1, j + 1, line[j] - 'A' + 1);
		}
	}
}

int main() {
	// freopen("E://input.txt", "r", stdin);
	int kase = 0;
	while (scanf("%s", line) != EOF) {
		++kase;
		if (kase != 1) {
			printf("\n");
		}
		memset(board, 0, sizeof board);
		memset(ans, 0, sizeof ans);
		process(0);
		for (int i = 1; i < 16; i++) {
			scanf("%s", line);
			process(i);
		}
		dlx.init(4096, 1024);
		for (int i = 1; i <= 4096; i++) {
			for (int j = 1; j <= 1024; j++) {
				if (board[i][j] == 1) {
					dlx.addlink(i, j);
				}
			}
		}

		dlx.dance(0);
		for (int i = 0; i < d; i++) {
			int num = (ans[i] - 1) / 16 + 1, val = ans[i] - 16 * (num - 1);
			int row = (num - 1) / 16 + 1, col = num - 16 * (row - 1);
			grid[row - 1][col - 1] = 'A' + val - 1;
		}
		for (int i = 0; i < 16; i++) {
			printf("%s\n", grid[i]);
		}
	}

	return 0;
}
```

# 20.2.18 HDU-4069 奇怪的数独 DFS+精确覆盖
给定的数独宫格并不是方方正正的了,而是奇形怪状的。先用dfs求一遍数独宫格号然后就套模板即可。

---
这题一开始TLE了一次,因此我把所有结果都搜了,其实只要搜到两种答案就不符合题意,即可返回;然后又WA了几次,这是因为没有再第一次搜索成功时就保存答案,而接下来的搜索就打断了原有的答案数组。

还有一个优化的地方,一开始1700ms,2.5M的内存,是因为我单独开了一个board数组用于存储Dancing Link的矩阵,其实不需要,只需要在insert的时候直接调用dlx.addlink插入就可以,因为这样也可以保证我插入的元素一定是竖着的最下面那个,横着的最右边那个。
```
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

const int fac[] = {128, 64, 32, 16};
const int go[][2] = {0,-1, 1,0, 0,1, -1,0};
const int maxM = 81 * 9 + 5;
const int maxN = 81 * 4 + 5;
const int maxMN = maxM * maxN;

int grid[10][10], cellNum[10][10];
int ans[maxM];

struct DLX {
	int m, n;
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxM], S[maxN];
	int sz;
	int solutions;

	void init(int m_, int n_) {
		m = m_; n = n_;
		for (int i = 0; i <= n; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = n, R[n] = 0;
		memset(H, -1, sizeof H);
		sz = n + 1;
		solutions = 0;
	}

	void addlink(int r, int c) {
		Row[sz] = r, Col[sz] = c;
		S[c]++;

		D[U[c]] = sz, U[sz] = U[c];
		U[c] = sz, D[sz] = c;

		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz, L[sz] = L[H[r]];
			R[sz] = H[r], L[H[r]] = sz;
		}
		sz++;
	}

	void remove(int c) {
		R[L[c]] = R[c], L[R[c]] = L[c];
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = U[j];
				D[U[j]] = D[j];
				S[Col[j]]--;
			}
		}
	}

	void resume(int c) {
		R[L[c]] = L[R[c]] = c;
		for (int i = D[c]; i != c; i = D[i]) {
			for (int j = R[i]; j != i; j = R[j]) {
				U[D[j]] = D[U[j]] = j;
				S[Col[j]]++;
			}
		}
	}

	void dance(int depth) {
		if (R[0] == 0) {
			solutions++;
			for (int i = 0; i < depth; i++) {
				int cell = (ans[i] - 1) / 9 + 1, val = ans[i] - 9 * (cell - 1);
				int row = (cell - 1) / 9 + 1, col = cell - 9 * (row - 1);
				grid[row - 1][col - 1] = val;
			}
			return;
		}

		int c = R[0];
		for (int i = c; i != 0; i = R[i]) {
			if (S[i] < S[c]) {
				c = i;
			}
		}

		remove(c);
		for (int i = D[c]; i != c; i = D[i]) {
			ans[depth] = Row[i];
			for (int j = R[i]; j != i; j = R[j]) {
				remove(Col[j]);
			}
			dance(depth+1);
			if (solutions >= 2) {
				return;
			}
			for (int j = L[i]; j != i; j = L[j]) {
				resume(Col[j]);
			}
		}
		resume(c);
	}
}dlx;

void dfs(int x, int y, int cell) {
	cellNum[x][y] = cell;
	int dir[4] = {0, 0, 0, 0};
	for (int i = 0; i < 4; i++) {
		if (grid[x][y] >= fac[i]) {
			grid[x][y] -= fac[i];
			dir[i] = 1;
		}
	}
	for (int i = 0; i < 4; i++) {
		if (dir[i]) continue;
		int nx = x + go[i][0];
		int ny = y + go[i][1];
		if (!cellNum[nx][ny]) {
			dfs(nx, ny, cell);
		}
	}
}

void insert(int x, int y, int val) {
	dlx.addlink(81 * x + 9 * y + val, 9 * x + y + 1);
	dlx.addlink(81 * x + 9 * y + val, 81 + 9 * x + val);
	dlx.addlink(81 * x + 9 * y + val, 81 * 2 + 9 * y + val);
	dlx.addlink(81 * x + 9 * y + val, 81 * 3 + 9 * (cellNum[x][y] - 1) + val);
}

int main() {
	// freopen("E://input.txt", "r", stdin);
	int t, kase = 0;
	scanf("%d", &t);
	while (t--) {
		kase++;
		memset(cellNum, 0, sizeof cellNum);
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				scanf("%d", &grid[i][j]);
			}
		}

		int cell = 1;
		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				if (!cellNum[i][j]) {
					dfs(i, j, cell++);
				}
			}
		}

		dlx.init(81 * 9, 81 * 4);

		for (int i = 0; i < 9; i++) {
			for (int j = 0; j < 9; j++) {
				if (grid[i][j] == 0) {
					for (int val = 1; val <= 9; val++) {
						insert(i, j, val);
					}
				}
				else {
					insert(i, j, grid[i][j]);
				}
			}
		}

		dlx.dance(0);

		printf("Case %d:\n", kase);
		if (dlx.solutions == 0) {
			printf("No solution\n");
		}
		else if (dlx.solutions > 1) {
			printf("Multiple Solutions\n");
		}
		else {
			for (int i = 0; i < 9; i++) {
				for (int j = 0; j < 9; j++) {
					printf("%d", grid[i][j]);
				}
				printf("\n");
			}
		}
	}
}
```



# 20.2.19 HDU 3335 Divisibility 重复覆盖解最大覆盖集 or 二分图匹配

给定n个数，当从中选择一个数a时，当能被a整除或者a整除的数都不能再选了，问最多能选多少个

---

重复覆盖求解,行和列都是n个数,当nums[i]和nums[j]存在整除关系时建边;重复覆盖是利用IDA*求取的最小搜索深度,而题目是最多选取,因此答案是重复覆盖解集合中的最大的那一个(而不是利用迭代加深搜索去寻找最小的那一个);因此剪枝暂时没用,也通过了。

---

重复覆盖的一个解集合S满足，其去掉任意一个S内的元素s，剩下的集合都不能满足覆盖的要求。

```
#include <cstdio>
#include <cstring>
#include <algorithm>

using namespace std;

const int maxM = 1005, maxN = 1005;
const int maxMN = maxM * maxN;

typedef long long ll;
ll nums[maxM];
int ans;

struct DLX {
	const int head = 0;
	int m, n;
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxM], S[maxN];
	int sz;

	void init(int m_, int n_) {
		m = m_, n = n_;
		for (int i = 0; i <= n; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = n, R[n] = 0;
		memset(H, -1, sizeof H);
		sz = n + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r, Col[sz] = c;
		S[c]++;

		D[U[c]] = sz, U[sz] = U[c];
		D[sz] = c, U[c] = sz;

		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz, L[sz] = L[H[r]];
			R[sz] = H[r], L[H[r]] = sz;
		}
		sz++;
	}

	void repeat_remove(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = L[i], R[L[i]] = R[i];
		}
	}

	void repeat_resume(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = R[L[i]] = i;
		}
	}

	int f() {
		bool vis[maxN];
		int res = 0;
		for (int i = R[head]; i != head; i = R[i]) vis[i] = true;
		for (int c = R[head]; c != head; c = R[c]) {
			if (vis[c]) {
				vis[c] = false;
				res++;
				for (int i = D[c]; i != c; i = D[i]) {
					for (int j = R[i]; j != i; j = R[j]) {
						vis[Col[j]] = false;
					}
				}
			}
		}
		return res;
	}

	void dance(int depth) {
		if (R[head] == head) {
			// depth是一个合法解
			ans = max(ans, depth);
			return;
		}
		// A*剪枝
		//if (f() + depth > ans) {
		//	return;
		//}

		int c = R[head];
		for (int i = R[head]; i != head; i = R[i]) {
			if (S[i] < S[c]) {
				c = i;
			}
		}

		for (int i = D[c]; i != c; i = D[i]) {
			repeat_remove(i);
			for (int j = R[i]; j != i; j = R[j]) {
				repeat_remove(j);
			}
			dance(depth+1);
			for (int j = L[i]; j != i; j = L[j]) {
				repeat_resume(j);
			}
			repeat_resume(i);
		}
	}
}dlx;

int main() {
	// freopen("E://input.txt", "r", stdin);
	int t;
	scanf("%d", &t);
	while (t--) {
		int n;
		scanf("%d", &n);
		for (int i = 1; i <= n; i++) {
			scanf("%lld", &nums[i]);
		}
		dlx.init(n, n);
		for (int i = 1; i <= n; i++) {
			for (int j = 1; j <= n; j++) {
				if (nums[i] % nums[j] == 0 || nums[j] % nums[i] == 0) {
					dlx.addlink(i, j);
				}
			}
		}
		ans = 0;
		dlx.dance(0);
		printf("%d\n", ans);
	}
	return 0;
}
```

# 20.2.20 HDU 4979 重复覆盖+线下打表+二进制集合枚举子集(★★★★)

从1到n中选取m个数和r个数,m>=r,问至少选多少个m个数的排列,可以覆盖到r个数的全排列;

---

典型重复覆盖问题,每一行是m个数的一个排列,每一列是r个数的一个排列,问题在于如何建表;转化为一个数学问题则是
```
从n个数中选取一个m个数的排列和一个r个数的排列,问r排列是否是m排列的子集
```
## 二进制集合枚举子集

这里用到了一个trick,即二进制集合枚举子集,我们用二进制来模拟每一个数的选择与否,那么对于一种选择状态s,其子集可以这样遍历
```
for (int i = s; i > 0; i = (i-1)&s)
```
每一次循环的i就是s的一个子集,这样理解:由于最后是&s因此s中为0的位在子集中也不可能为1,因此不考虑0,i-1就是将原来是1的位不断-1。

---

回到该题上,我们首先保存状态i的1的个数cnt(代表选了几个数),然后若cnt为m,则按照上面求取子集的方式建边,这里又涉及到一个问题,由于i是从大到小遍历的,因此在DLX建边中水平方向要插到第一个节点的后面而不是最后一个节点的后面;

最后就是这题需要离线打大表,无论怎么剪枝搜索都是时间爆掉,本地将全部答案求出来,打表,然后直接访问输出

还有一个点就是IDA*求重复覆盖时dance函数的写法;我一开始是这样写的
```
for (int maxd = 0; ; maxd++) {
    if (dlx.dance(0, maxd)) {
        break;
    }
}
```
这种就是从0开始一直搜一直搜直到找到maxd答案,但是今天我看到了是这样写的,即一开始设置最大深度为无穷大,一旦找到一个解,若该解比最大深度小则更新最大深度;剪枝的判断就是当depth+A() >= Min即可,意思就是当前情况往下接着搜也找不到比我目前发现的Min更小的答案了,因此你就别搜了;就本题而言,第二种方法更快;

```
#include <cstdio>
#include <cstring>

const int MAXN = 8;
const int MAXR = 1000;
const int MAXC = 1000;
const int MAXNODE = MAXR * MAXC;
const int INF = 0x3f3f3f3f;

int stateInCol[MAXC], ccnt;
int rcnt;
int bitcnt[1 << MAXN];
int C[MAXN + 1][MAXN + 1];

struct DLX
{
	int sz;
	int H[MAXR], S[MAXC];
	int row[MAXNODE], col[MAXNODE];
	int U[MAXNODE], D[MAXNODE], L[MAXNODE], R[MAXNODE];
	int Min;

	void Init(int n)
	{
		for (int i = 0; i <= n; ++i)
		{
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
		}
		L[0] = n;
		R[n] = 0;

		sz = n + 1;
		memset(S, 0, sizeof(S));
		memset(H, -1, sizeof(H));
	}

	void Link(const int& r, const int& c)
	{
		row[sz] = r;
		col[sz] = c;
		D[sz] = D[c];
		U[D[c]] = sz;
		D[c] = sz;
		U[sz] = c;
		if (H[r] == -1)
		{
			H[r] = L[sz] = R[sz] = sz;
		}
		else
		{
		    // 注意这个地方是插入第一个结点的后面
		    // 由于第一个节点是自己指向自己,因此这种插法就是
		    // 尾插法
			R[sz] = R[H[r]];
			L[R[H[r]]] = sz;
			R[H[r]] = sz;
			L[sz] = H[r];
		}
		S[c]++;
		sz++;
	}

	void Remove(const int& c)
	{
		for (int i = D[c]; i != c; i = D[i])
		{
			L[R[i]] = L[i];
			R[L[i]] = R[i];
		}
	}

	void Restore(const int& c)
	{
		for (int i = U[c]; i != c; i = U[i])
		{
			L[R[i]] = i;
			R[L[i]] = i;
		}
	}

	int A()
	{
		int ret = 0;
		bool vis[MAXC];

		memset(vis, 0, sizeof(vis));
		for (int i = R[0]; i != 0; i = R[i])
		{
			if (!vis[i])
			{
				vis[i] = true;
				++ret;
				for (int j = D[i]; j != i; j = D[j])
				{
					for (int k = R[j]; k != j; k = R[k])
					{
						vis[col[k]] = true;
					}
				}
			}
		}

		return ret;
	}

	void Dfs(int cur)
	{
	    // 留意DLX这种IDA*的写法
		if (cur + A() >= Min) return;

		if (R[0] == 0)
		{
			if (cur < Min)
			{
				Min = cur;
			}
			return;
		}

		int c = R[0];
		for (int i = R[0]; i != 0; i = R[i])
		{
			if (S[i] < S[c])
			{
				c = i;
			}
		}

		for (int i = D[c]; i != c; i = D[i])
		{
			Remove(i);
			for (int j = R[i]; j != i; j = R[j])
			{
				Remove(j);
			}
			Dfs(cur + 1);
			for (int j = L[i]; j != i; j = L[j])
			{
				Restore(j);
			}
			Restore(i);
		}
	}

	int Solve()
	{
		Min = INF;
		Dfs(0);
		return Min;
	}

} dlx;
// LeetCode 求解数x二进制上有多少1 x &= (x-1)
int Bitcnt(int x)
{
	int ret = 0;

	while (x)
	{
		ret += (x & 1);
		x >>= 1;
	}

	return ret;
}

void GetBitcnt()
{
	for (int i = 0; i < (1 << MAXN); ++i)
	{
		bitcnt[i] = Bitcnt(i);
	}
}
// 计算组合数Cnm
void GetC()
{
	for (int i = 1; i <= MAXN; ++i)
	{
		C[i][0] = C[i][i] = 1;
		for (int j = 1; j < i; ++j)
		{
			C[i][j] = C[i - 1][j] + C[i - 1][j - 1];
		}
	}
}

void Init()
{
	GetBitcnt();
	GetC();
}

void Solve(int N, int M, int R)
{
	ccnt = 0;
	// 首先将r个数的排列标列号,从小到大
	for (int i = 1; i < (1 << N); ++i)
	{
		if (bitcnt[i] == R)
		{
			stateInCol[i] = ++ccnt;
		}
	}

	dlx.Init(C[N][R]);

	rcnt = 0;
	for (int i = 1; i < (1 << N); ++i)
	{
		if (bitcnt[i] == M)
		{
			++rcnt;
			// 遍历二进制j的全部子集
			for (int j = i; j > 0; j = (i & (j - 1)))
			{
				if (bitcnt[j] == R)
				{
					dlx.Link(rcnt, stateInCol[j]);
				}
			}
		}
	}

	printf("%d", dlx.Solve());
}

// 离线打表,本地运行打印出来后直接赋给数组
void SaveTable()
{
	puts("{");
	for (int N = 1; N <= MAXN; ++N)
	{
		puts("  {");
		for (int M = 1; M <= N; ++M)
		{
			printf("    {");
			for (int R = 1; R <= M; ++R)
			{
				if (R > 1)
				{
					printf(", ");
				}
				Solve(N, M, R);
			}
			printf("}");
			if (M == N)
			{
				puts("");
			}
			else
			{
				puts(",");
			}
		}
		printf("  }");
		if (N == MAXN)
		{
			puts("");
		}
		else
		{
			puts(",");
		}
	}
	puts("}");
}
// 这是solveTable的结果,自己电脑上运行了好几分钟
int st[8][8][8] = {
	{
		{ 1 }
	},
	{
		{ 2 },
		{ 1, 1 }
	},
	{
		{ 3 },
		{ 2, 3 },
		{ 1, 1, 1 }
	},
	{
		{ 4 },
		{ 2, 6 },
		{ 2, 3, 4 },
		{ 1, 1, 1, 1 }
	},
	{
		{ 5 },
		{ 3, 10 },
		{ 2, 4, 10 },
		{ 2, 3, 4, 5 },
		{ 1, 1, 1, 1, 1 }
	},
	{
		{ 6 },
		{ 3, 15 },
		{ 2, 6, 20 },
		{ 2, 3, 6, 15 },
		{ 2, 3, 4, 5, 6 },
		{ 1, 1, 1, 1, 1, 1 }
	},
	{
		{ 7 },
		{ 4, 21 },
		{ 3, 7, 35 },
		{ 2, 5, 12, 35 },
		{ 2, 3, 5, 9, 21 },
		{ 2, 3, 4, 5, 6, 7 },
		{ 1, 1, 1, 1, 1, 1, 1 }
	},
	{
		{ 8 },
		{ 4, 28 },
		{ 3, 11, 56 },
		{ 2, 6, 14, 70 },
		{ 2, 4, 8, 20, 56 },
		{ 2, 3, 4, 7, 12, 28 },
		{ 2, 3, 4, 5, 6, 7, 8 },
		{ 1, 1, 1, 1, 1, 1, 1, 1 }
	}
};
// 直接打表AC
int main()
{
	int t, kase = 0;
	scanf("%d", &t);
	while (t--) {
		int n, m, r;
		scanf("%d %d %d", &n, &m, &r);
		printf("Case #%d: %d\n", ++kase, st[n-1][m-1][r-1]);
	}

	return 0;
}
```
PS：龟龟,原来这才是真正的“打表”,我之前以为的打表都弱爆了!

# 20.2.21 HDU 5046 Airport 最大值最小值二分+重复覆盖

给定n个城市，从中选取k个城市作机场，令di代表第i个城市到最近机场的距离，现在求max(di, i=1,..,n)的最小值

## 最大值最小化问题

首先是一个最大值最小化的问题,这类问题一般都用二分法解决,比如
```
把一个包含n个数的序列划分为连续的m个子序列,设第i个序列的各数之和为S(i),求所有S(i)的最大值最小是多少？
```
对于一次划分,求一个x,使得对任意的S(i),都有S(i) <= x;这个条件保证了x是所有S(i)中的最大值;我们需要求的就是满足该条件的最小的x;
则思路为二分搜索x,对于当前的值x,判断是否可以找到一次划分满足S(i) <= x,如果满足,真实答案一定不比x要大,收缩右边界；如果不满足,说明找不到任何一种划分,使其序列和的最大值小于x,那么更不用提比当前x还小的值;因此二分调性成立；

回到这道题,进行同样的转化,我们也是寻找一个距离x,该距离满足max{di, i=1..n},现在的目标就是搜索最小的x；
思路同样二分搜索x,对于当前的值x,判断是否可以找到一次选取满足所有的di <= x,即对于任意城市到最近机场的距离都小于x；
转化为重复覆盖问题,行列都为城市,当i,j城市之间的距离<=x时建边，问搜索深度为k时是否有解;

```C++
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <cmath>

using namespace std;

const int maxM = 100, maxN = 100;
const int maxMN = maxM * maxN;

struct DLX {
	const int head = 0;
	int m, n;
	int U[maxMN], D[maxMN], L[maxMN], R[maxMN];
	int Row[maxMN], Col[maxMN];
	int H[maxM], S[maxN];
	int sz;

	void init(int m_, int n_) {
		m = m_, n = n_;
		for (int i = 0; i <= n; i++) {
			U[i] = D[i] = i;
			L[i] = i - 1;
			R[i] = i + 1;
			S[i] = 0;
		}
		L[0] = n, R[n] = 0;
		memset(H, -1, sizeof H);
		sz = n + 1;
	}

	void addlink(int r, int c) {
		Row[sz] = r, Col[sz] = c;
		S[c]++;

		D[U[c]] = sz, U[sz] = U[c];
		U[c] = sz, D[sz] = c;

		if (H[r] < 0) {
			H[r] = L[sz] = R[sz] = sz;
		}
		else {
			R[L[H[r]]] = sz, L[sz] = L[H[r]];
			L[H[r]] = sz, R[sz] = H[r];
		}
		sz++;
	}

	void repeat_remove(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = L[i];
			R[L[i]] = R[i];
		}
	}

	void repeat_resume(int c) {
		for (int i = D[c]; i != c; i = D[i]) {
			L[R[i]] = R[L[i]] = i;
		}
	}

	int f() {
		int res = 0;
		bool vis[maxN];
		memset(vis, false, sizeof vis);
		for (int c = R[head]; c != head; c = R[c]) {
			if (!vis[c]) {
				vis[c] = true;
				res++;
				for (int i = D[c]; i != c; i = D[i]) {
					for (int j = R[i]; j != i; j = R[j]) {
						vis[Col[j]] = true;
					}
				}
			}
		}
		return res;
	}

	bool dance(int depth, int maxd) {
		if (depth + f() > maxd) return false;
		if (R[head] == head) return true;
		int c = R[head];
		for (int i = c; i != head; i = R[i]) {
			if (S[i] < S[c]) {
				c = i;
			}
		}

		for (int i = D[c]; i != c; i = D[i]) {
			repeat_remove(i);
			for (int j = R[i]; j != i; j = R[j]) {
				repeat_remove(j);
			}
			if (dance(depth + 1, maxd)) return true;
			for (int j = L[i]; j != i; j = L[j]) {
				repeat_resume(j);
			}
			repeat_resume(i);
		}
		return false;
	}
}dlx;

typedef long long ll;
ll x[maxM], y[maxM];

ll getDistance(int i, int j) {
	return abs(x[i]-x[j]) + abs(y[i]-y[j]);
}
int main() {
	// freopen("E://input.txt", "r", stdin);
	int t, kase = 0;
	scanf("%d", &t);
	while (t--) {
		kase++;
		int n, k;
		scanf("%d %d", &n, &k);
		for (int i = 1; i <= n; i++) {
			scanf("%lld %lld", &x[i], &y[i]);
		}
		
		ll l = 0, h = 1e10, mid = 0;
		while (l < h) {
			mid = (l + h) / 2;
			dlx.init(n, n);
			for (int i = 1; i <= n; i++) {
				for (int j = 1; j <= n; j++) {
					if (getDistance(i, j) <= mid) {
						dlx.addlink(i, j);
					}
				}
			}
			int res = dlx.dance(0, k);
			if (res) {
				// 答案肯定不会比mid大
				h = mid;
			}
			else {
				// 答案一定比mid大
				l = mid+1;
			}
		}
		printf("Case #%d: %lld\n", kase, l);
	}
	return 0;
}
```