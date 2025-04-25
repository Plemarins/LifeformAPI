PROGRAM ArtificialLife
  IMPLICIT NONE
  INTEGER, PARAMETER :: N = 50        ! グリッドサイズ
  INTEGER, PARAMETER :: TMAX = 100    ! 最大時間ステップ
  INTEGER, PARAMETER :: NPOP = 20     ! GAの個体数
  REAL, PARAMETER :: D = 0.1          ! 拡散係数
  REAL, PARAMETER :: ALPHA = 0.05     ! 減衰率
  REAL, PARAMETER :: MUT_RATE = 0.01  ! 突然変異率
  INTEGER :: i, j, t, k
  REAL :: u(N,N), u_new(N,N)          ! セルの状態
  REAL :: w(9,NPOP)                   ! NNの重み（9近傍）
  REAL :: fitness(NPOP)               ! 適応度
  REAL :: u_enemy(N,N)                ! 敵対的セルの状態（簡略化）

  ! 初期化
  CALL RANDOM_SEED()
  DO i = 1, N
    DO j = 1, N
      CALL RANDOM_NUMBER(u(i,j))
      u_enemy(i,j) = u(i,j) * 0.5  ! 敵対的セルの初期化
    END DO
  END DO
  DO k = 1, NPOP
    DO i = 1, 9
      CALL RANDOM_NUMBER(w(i,k))
      w(i,k) = w(i,k) - 0.5  ! [-0.5, 0.5] の範囲
    END DO
  END DO

  ! メインループ
  DO t = 1, TMAX
    ! CA更新（反応拡散 + NN）
    DO i = 2, N-1
      DO j = 2, N-1
        ! 近傍状態の取得
        REAL :: nb(9)
        nb(1) = u(i-1,j-1); nb(2) = u(i-1,j); nb(3) = u(i-1,j+1)
        nb(4) = u(i,j-1);   nb(5) = u(i,j);   nb(6) = u(i,j+1)
        nb(7) = u(i+1,j-1); nb(8) = u(i+1,j); nb(9) = u(i+1,j+1)
        
        ! NN出力（代表個体の重みを使用）
        REAL :: nn_out
        nn_out = 0.0
        DO k = 1, 9
          nn_out = nn_out + w(k,1) * nb(k)
        END DO
        nn_out = 1.0 / (1.0 + EXP(-nn_out))  ! シグモイド
        
        ! 拡散項
        REAL :: laplacian
        laplacian = u(i+1,j) + u(i-1,j) + u(i,j+1) + u(i,j-1) - 4.0 * u(i,j)
        
        ! ミニマックス簡略化（敵対的セルの影響）
        REAL :: mm_term
        mm_term = MIN(u(i,j), u_enemy(i,j))  ! 簡略化された報酬
        
        ! 状態更新
        u_new(i,j) = u(i,j) + D * laplacian + nn_out - ALPHA * u(i,j) + 0.1 * mm_term
        IF (u_new(i,j) < 0.0) u_new(i,j) = 0.0
        IF (u_new(i,j) > 1.0) u_new(i,j) = 1.0
      END DO
    END DO
    u = u_new

    ! GAによる進化（簡略化）
    DO k = 1, NPOP
      ! 適応度計算（例：グリッド全体の生存率）
      fitness(k) = SUM(u) / (N * N)
      ! ノイズを追加して個体差を模擬
      CALL RANDOM_NUMBER(nn_out)
      fitness(k) = fitness(k) + 0.01 * nn_out
    END DO
    
    ! 選択と突然変異
    CALL SelectionAndMutation(w, fitness, NPOP, MUT_RATE)

    ! 結果出力（簡略化）
    PRINT *, 'Time = ', t, ' Avg State = ', SUM(u) / (N * N)
  END DO

CONTAINS
  SUBROUTINE SelectionAndMutation(w, fitness, npop, mut_rate)
    REAL, INTENT(INOUT) :: w(9, npop)
    REAL, INTENT(IN) :: fitness(npop), mut_rate
    INTEGER :: i, k
    REAL :: w_new(9, npop)
    ! 簡略化：適応度の高い個体をコピー
    DO k = 1, npop
      w_new(:,k) = w(:,1)  ! 最良個体をコピー
      DO i = 1, 9
        CALL RANDOM_NUMBER(nn_out)
        IF (nn_out < mut_rate) THEN
          CALL RANDOM_NUMBER(nn_out)
          w_new(i,k) = w_new(i,k) + 0.1 * (nn_out - 0.5)
        END IF
      END DO
    END DO
    w = w_new
  END SUBROUTINE
END PROGRAM
