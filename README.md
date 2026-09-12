# 证明不定方程 $2^n - n = k^2$ 的正整数解仅有 $n=1,7$

证明见 [proof.tex](proof.tex), [proof.pdf](proof.pdf), 解个数的有限性证明由 gpt-5.6-sol xhigh 生成, 检验 $n\le 5\times 10^5$ 情形由 GLM-5.3 max 生成，检验代码见 [check.py](check.py), 多进程版本见 [check_multiprocess.py](check_multiprocess.py)

扩展情形 $u^n+f(n)=k^2,u\geq 2, \log|f(n)|=o(n/\log n)$ 解个数的有限性证明见 [proof_finiteness_u_n_plus_f_n.pdf](proof_finiteness_u_n_plus_f_n.pdf)，它由 GLM-5.3 在上述证明基础上生成
