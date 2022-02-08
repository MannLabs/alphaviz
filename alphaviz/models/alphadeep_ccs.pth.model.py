import torch
import alphadeep.model.base as model_base
class Model(torch.nn.Module):
    def __init__(self,
        dropout=0.1
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        hidden = 256

        self.ccs_encoder = \
            model_base.Input_AA_CNN_LSTM_cat_Charge_Encoder(
                hidden
            )

        self.ccs_decoder = model_base.LinearDecoder(
            hidden+1, 1
        )

    def forward(self,
        aa_indices,
        mod_x,
        charges,
    ):
        x = self.ccs_encoder(aa_indices, mod_x, charges)
        x = self.dropout(x)
        x = torch.cat((x, charges),1)
        return self.ccs_decoder(x).squeeze(1)
