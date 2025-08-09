use dqmodule
implicit none
type(dq_real)    :: x
character(128)   :: s
x = dqpi()                          ! 円周率
call dqwrite(6, 80, 65, x)          ! stdout に E80.65 で出力（幅は 65+10 以上）  [README 記載]
s = '3.14159265358979323846264338327950288419716939937510'
x = dqreal(s)                       ! 文字列から厳密に DQR へ
call dqwrite(6, 40, 30, x)
end
