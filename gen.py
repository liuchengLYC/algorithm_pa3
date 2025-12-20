#!/usr/bin/env python3
import random
from pathlib import Path

# ============================
# Configuration
# ============================
xSize = 1000
ySize = 1000
num_layers = 2
wlViaCost = 100
num_nets = 2000

cap_output = "case7.cap"
net_output = "case7.net"

rng = random.Random(12345)


# ============================
# Generate .cap file
# ============================
def generate_cap():
    with open(cap_output, "w") as f:
        # 第一行：<layers> <xSize> <ySize>
        f.write(f"{num_layers} {xSize} {ySize}\n")

        # 第二行：wlViaCost
        f.write(f"{wlViaCost}\n")

        # 第三行：horizontal distance（長度 xSize-1）
        f.write(" ".join(["6000"] * (xSize - 1)) + "\n")

        # 第四行：vertical distance（長度 ySize-1）
        f.write(" ".join(["5700"] * (ySize - 1)) + "\n")

        # Layer 0
        f.write("Metal1 H\n")
        for _ in range(ySize):
            row = [str(rng.randint(1, 2)) for _ in range(xSize)]
            f.write(" ".join(row) + "\n")

        # Layer 1
        f.write("Metal2 V\n")
        for _ in range(ySize):
            row = [str(rng.randint(1, 2)) for _ in range(xSize)]
            f.write(" ".join(row) + "\n")

    print(f"[OK] 生成 {cap_output}")


# ============================
# Generate .net file
# ============================
def generate_net():
    with open(net_output, "w") as f:
        for i in range(num_nets):
            f.write(f"net{i}\n")
            f.write("(\n")

            # 隨機 pin1
            l1 = rng.randint(0, 1)
            c1 = rng.randrange(xSize)
            r1 = rng.randrange(ySize)

            # 隨機 pin2（確保不同位置）
            while True:
                l2 = rng.randint(0, 1)
                c2 = rng.randrange(xSize)
                r2 = rng.randrange(ySize)
                if (l2, c2, r2) != (l1, c1, r1):
                    break

            f.write(f"({l1},{c1},{r1})\n")
            f.write(f"({l2},{c2},{r2})\n")
            f.write(")\n")

    print(f"[OK] 生成 {net_output}")


# ============================
# Main
# ============================
if __name__ == "__main__":
    generate_cap()
    generate_net()
    print("[DONE] case7 資料生成完成！")
